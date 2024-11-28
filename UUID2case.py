import os

manifest_file = "gdc_manifest.2024-11-15.txt"

# join strings of ids
ids = [" "]
manifest_length = 0
with open(manifest_file, "r", encoding="utf-8") as f:
    fields = f.readline()
    while True:
        if fields:
            ids.append("\"{}\"".format(fields.strip().split("\t")[0]))
        fields = f.readline()
        manifest_length += 1
        
        if not fields:
            break
        ids.append(", ")
file_ids = "".join(ids)

# prepare Payload.txt
Part1= '{\"filters\":{\"op\":\"in\",\"content\":{\"field\":\"files.file_id\",\"value\":[ '
Part2= '] }},\"format\":\"TSV\",\"fields\":\"cases.project.project_id,cases.samples.tissue_type,file_id\",\"size\":'
Part3= "".join(["\"", str(manifest_length), "\"", "}"])
Sentence = " ".join([Part1, file_ids, Part2, Part3])

with open("Payload.txt", "w", encoding="utf-8") as f:
    f.write(Sentence)

# cmd execute
cmd = "curl --request POST --header \"Content-Type: application/json\" --data @Payload.txt \"https://api.gdc.cancer.gov/files\" > File_metadata.txt"
os.system(cmd)