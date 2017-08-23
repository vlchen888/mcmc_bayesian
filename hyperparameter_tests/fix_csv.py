import csv
def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

methods = ['tau_alpha', 'tau_M', 'alpha_M', 'tau_alpha_M']

with open('record.csv','r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    with open('new_record.csv','w') as newfile:
        writer = csv.writer(newfile, delimiter=',')
        trial = 0
        method = ''
        for row in reader:
            if is_number(row[0]):
                trial = int(row[0])
                method = methods.index(row[1])+1
            elif row[0] is '':
                sample = row[2]
                accuracy = row[3]
                writer.writerow([trial] + [method] + [sample] + [accuracy])