import json

count = 0
project = dict()
with open("File_metadata.txt", "r", encoding="utf-8") as file:
    line = file.readline()
    while True:
        line = file.readline()
        if not line:
            print("number of entities:", count)
            break
        count += 1
        fields = line.strip().split("\t")
        if fields[1] == "Normal":
            project[fields[2]] = "Normal"
        elif fields[1] == "Tumor":
            tcgaid = fields[0].split("-") 
            project[fields[2]] = tcgaid[1]
        else:
            raise ValueError("unknown tissue type:", fields[1])
with open("project_map.json", "w", encoding="utf-8") as f:
    json.dump(project, f, ensure_ascii=False, indent=4)
print("save as .json file")