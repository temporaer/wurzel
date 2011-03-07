#include <opencv/cv.h>
//#include <opencv/highgui.h>
#include "cvhough.cpp"
#include <math.h>

float*
find_circles(const std::vector<point_type>& rec, unsigned int resx, unsigned int resy, const stat_t& dim0, const stat_t& dim1, const stat_t& dim2)
{
    IplImage* img = cvCreateImage(cvSize(resx,resy), IPL_DEPTH_32F, 1);

    IplImage* gray = cvCreateImage(cvSize(resx,resy), 8, 1);
    cvSetZero(gray);
    int maxv = -1;
    foreach(const point_type& p, rec){
	    // scaling is according to dim0 only!
	   unsigned int px = (unsigned int)(resx * (p[0]-acc::min(dim0))/(acc::max(dim0)-acc::min(dim0)));
	   unsigned int py = (unsigned int)(resy * (p[1]-acc::min(dim1))/(acc::max(dim0)-acc::min(dim0)));
	   float v = CV_IMAGE_ELEM(img, float, px, py);

	   // ignore floor
	   if(p[2]<acc::min(dim2) + 0.1 * (acc::max(dim2)-acc::min(dim2)))
		   continue;
	   //CV_IMAGE_ELEM(gray, unsigned char, px, py) = (v==250) ? 250 :(v+50);
	   //CV_IMAGE_ELEM(gray, unsigned char, px, py) = 255;
	   CV_IMAGE_ELEM(img, float, px, py) +=1;
	   maxv = std::max(maxv, (int)v+1);
    }
    std::cout << "Maximum v: "<<maxv<<std::endl;
    std::ofstream os_img("gray.txt");
    for (int i = 0; i < resx; ++i) {
	    for (int j = 0; j < resy; ++j) {
		    float v = CV_IMAGE_ELEM(img, float, i, j);
		    CV_IMAGE_ELEM(gray, unsigned char, i, j) = (unsigned char) (v/maxv*255);
		    os_img<<(int)(v/maxv*255)<<" ";
	    }
	    os_img<<std::endl;
    }
    CvMemStorage* storage = cvCreateMemStorage(0);
    cvSmooth(gray, gray, CV_GAUSSIAN, 9, 9); 
    CvSeq* circles = cvHoughCircles_NoCanny(gray, storage, 
        CV_HOUGH_GRADIENT, 2, gray->height/100, 200, 100,0,0);

    std::cout << "Detected "<<circles->total<<" circles"<<std::endl;
    double max_radius=-1;
    int arg_max_radius=0;
    for (int i = 0; i < circles->total; i++) 
    {
         float* p = (float*)cvGetSeqElem( circles, i );
	 //cvCircle( gray, cvPoint(cvRound(p[0]),cvRound(p[1])), 0, 255, -1, 8, 0 );
	 //cvCircle( gray, cvPoint(cvRound(p[0]),cvRound(p[1])), cvRound(p[2]), 255, 3, 8, 0 );
	 //cvCircle( img, cvPoint(cvRound(p[0]),cvRound(p[1])), 3, CV_RGB(0,255,0), -1, 8, 0 );
	 //cvCircle( img, cvPoint(cvRound(p[0]),cvRound(p[1])), cvRound(p[2]), CV_RGB(255,0,0), 3, 8, 0 );

	 // make sure circle is completely inside image
	 if(p[0]<p[2]) continue;
	 if(p[1]<p[2]) continue;

	 if(p[0]>resx-p[2]) continue;
	 if(p[1]>resy-p[2]) continue;
	 std::cout << "Circle: "<<cvRound(p[0])<<", "<<cvRound(p[1])<<";  "<< cvRound(p[2])<<std::endl;
	 if(p[2]>max_radius) {
		 arg_max_radius = i;
		 max_radius = p[2];
	 }
    }
    if(max_radius<0)
	    return NULL;
    float* best = (float*)cvGetSeqElem(circles, arg_max_radius);
    // backtransform
    best[0] /= resx;
    best[0] *= (acc::max(dim0)-acc::min(dim0));
    best[0] += acc::min(dim0);

    best[1] /= resy;
    best[1] *= (acc::max(dim0)-acc::min(dim0));
    best[1] += acc::min(dim1);

    best[2] /= resx;
    best[2] *= (acc::max(dim0)-acc::min(dim0));
    return best;
}
