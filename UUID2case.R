# prepare UUID
manifest= "gdc_manifest.2024-11-15.txt" #Manifest name 
x=read.table(manifest,header = T)
manifest_length= nrow(x)
id= toString(sprintf('"%s"', x$id))

# prepare Payload.txt
Part1= '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '
#Part2= '] }},"format":"TSV","fields":"file_id,file_name,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.project.project_id","size":'
Part2= '] }},"format":"TSV","fields":"cases.project.project_id,cases.samples.tissue_type,file_id","size":'
Part3= paste("\"",manifest_length, "\"", "}")
Sentence= paste(Part1,id,Part2,Part3, collapse=" ")
write.table(Sentence,"Payload.txt",quote=F,col.names=F,row.names=F)

# cmd execute
#curl --request POST --header "Content-Type: application/json" --data @Payload.txt "https://api.gdc.cancer.gov/files" > File_metadata.txt
