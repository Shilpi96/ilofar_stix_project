; Search AIA 171 Angstrom data for the 6th december 2010 with a cadence of 3600 seconds (1 hour)
results = VSO_SEARCH('6-dec-2010', '7-dec-2010', inst='aia', wave='171', sample=3600)

; Sort the results by date and print the date
results = results(SORT(cat.time.start))
print, results.time.start

; Download the first two from ROB
status = VSO_GET(results[0:2], site='rob', pixels=4096, filenames=files, /rice,  /use_network)

; The local filenames will be in the files variable
print, files
