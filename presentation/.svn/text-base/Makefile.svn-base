SOURCE=$(wildcard *.tex)
TARGET=$(patsubst %.tex,%,$(SOURCE))
INCLUDEONLYFRAMES=
WEBWIDTH=267
WEBHEIGHT=200
.PHONY : all clean distclean figures current handout dirs

all: dirs screen 

screen: dirs
	@perl -pi'.orig' -e 's/,compress,handout\]/,compress]/' $(SOURCE);
	@rubber -W all -s --into .mk -I pic -I package -I esann2010_lateral --pdf $(TARGET)
	@mv .mk/$(TARGET).pdf $(TARGET)_screen.pdf
handout: dirs
	perl -pi'.orig' -e 's/,compress\]/,compress,handout]/' $(SOURCE);
	rubber -W all -s --into .mk -I pic  -I package --pdf $(TARGET)
	mv .mk/$(TARGET).pdf $(TARGET)_handout.pdf
	#pdfnup --paper a4paper --orient landscape --nup 3x2 --frame true --outfile bla.pdf $(TARGET)_handout.pdf
	#mv bla.pdf $(TARGET)_handout.pdf
current: dirs
	perl -pi'.orig' -e 's/,compress,handout\]/,compress]/' $(SOURCE)
	perl -pi'.orig' -e 's/^% \\includeonly/\\includeonly/' $(SOURCE)
	rubber -W all -s -I pic --into .mk --pdf $(TARGET)
	mv .mk/$(TARGET).pdf $(TARGET)_current.pdf
	mv $(SOURCE).orig $(SOURCE)

web: handout screen
	acroread -toPostScript -expand -size $(WEBWIDTH)x$(WEBHEIGHT) $(TARGET)_screen.pdf
	acroread -toPostScript -shrink -landscape -size a4 $(TARGET)_handout.pdf
	ps2pdf14 -dPDFSETTINGS=/screen -dAutoFilterxxxImages=false -dDEVICEWIDTHPOINTS=$(WEBWIDTH) -dDEVICEHEIGHTPOINTS=$(WEBHEIGHT) -dFIXEDMEDIA $(TARGET)_screen.ps $(TARGET)_web.pdf
	ps2pdf14 -dPDFSETTINGS=/screen $(TARGET)_handout.ps $(TARGET)_handout_web.pdf

impressive: dirs
	#impressive-gettransitions $(TARGET).tex
	impressive --transition None -I presentation.pdf.info -g 1024x768 -e -d 25m -M -u 10 $(TARGET)_screen.pdf
impressive2: dirs
	impressive --cachefile .kjn/cache -c persistent -e -d 25m -M -u 10 --transition None $(TARGET)_screen.pdf
impressive_prepare: screen
	#impressive-gettransitions $(TARGET).tex
	impressive  --transition None -g 1024x768 -I presentation.pdf.info -d 25m -M -u 10 --cachefile .kjn/cache $(TARGET)_screen.pdf

%.eps : %.fig
	fig2dev -L eps $< $@

dirs:
	@mkdir -p .mk 
	@mkdir -p .kjn
	@cp package/* .mk/

clean:
	rm -rf .mk
distclean:clean
	rm -f $(TARGET).pdf
	rm -f $(TARGET)_*.pdf
	rm -f $(TARGET)_*.ps
