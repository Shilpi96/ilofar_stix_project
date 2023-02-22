; Download aia files with ssw/idl
results = VSO_SEARCH('11-nov-2022 11:40', '11-nov-2022 11:45', inst='aia', wave='171')

; Sort the results by date and print the date
results = results(SORT(cat.time.start))
print, results.time.start

; Download the first two from ROB
status = VSO_GET(results[0:2], site='rob', pixels=4096, filenames=files, /rice,  /use_network)

; The local filenames will be in the files variable
print, files
