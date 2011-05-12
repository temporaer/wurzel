import xml.etree.ElementTree as et

class config:
	def __init__(self,fn):
		dom       = et.parse(fn)
		root      = dom.getroot()
		datafiles = root.find("datafiles")
		D         = datafiles.findall("datafile")
		self.all_data = D
	def read_dtype(self, id):
		return self.get(id).get("read-dtype")
	def has_rohr(self, id):
		return self.get(id).find("has-rohr").text in ["true", "1"]
	def read_shape(self, id):
		return [int(self.get(id).find("read-shape").get(a)) for a in ["x", "y", "z"]]
	def shape(self, id):
		return [int(self.get(id).find("shape").get(a)) for a in ["x", "y", "z"]]
	def scale(self, id):
		return float(self.get(id).find("scale").text)
	def get(self, id):
		s = "../%s" % id
		for x in self.all_data:
			#print "compare: `%s'  with `%s'" %(id, x.find("base-name").text)
			if x.find("base-name").text == s:
				#print "success."
				return x
	def all_bases(self):
		return [x.find("base-name").text for x in self.all_data]

