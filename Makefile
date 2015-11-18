figure-labels.png: figure-labels.R data.by.type.RData
	R --no-save < $<
data.by.type.RData: data.by.type.R pairs.txt
	R --no-save < $<