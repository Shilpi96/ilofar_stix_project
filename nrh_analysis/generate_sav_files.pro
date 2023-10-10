;Go to swwidl and then type .r read_nrh_maps and then .r generate_sav_files and then again type generate_sav_files.
pro generate_sav_files
  
root_dir ="/mnt/LOFAR-PSP/ilofar_STIX_shilpi/all_radio_data/nrh_fits/"
cd, root_dir
files = findfile("*nrh*.fts")

for i=0, n_elements(files) do read_nrh_maps, files[i]
;read_nrh_maps,files[7]
end
