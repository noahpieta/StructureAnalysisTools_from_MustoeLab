import requests
import argparse
import xml.etree.ElementTree as ET
import os

# Function download eCLIP bed files from ENCODE database
def getdata(assembly, outpath=None):
    
    if os.path.isdir(outpath):
        print(f'{outpath} folder already exists. Files will be downloaded the same folder.')
    else:
        os.makedirs(outpath)
        print(f'Created {outpath} folder.')

    print("Searching ENCODE database for eCLIP experiments")
    accession = []
    bed_accession = []   

    # Force return from the server in JSON format
    headers = {'accept': 'application/json'}
    # GET accession number for experiments
    url_bed = 'https://www.encodeproject.org/search/?type=Experiment&assay_title=eCLIP' \
                  '&frame=object&limit=all&files.file_type=bed+narrowPeak'
    response = requests.get(url_bed, headers=headers)
    # Extract the json response as a Python dictionary and append accession to list
    experiment_dict = response.json()
    for entry in experiment_dict['@graph']:
        accession.append(entry['accession'])
       
    # GET accession number for IDR bed files, store in a dict.
    # Key is experiment acc and values are bed files acc for that experiment.
    for exp in accession:
        url2 = f'https://www.encodeproject.org/search/?type=File&dataset=/experiments/{exp}/' \
                   '&file_format=bed&format=json&frame=object&limit=all'
        response2 = requests.get(url2, headers=headers)
        bed_dict = response2.json()
        for i in bed_dict['@graph']:
            # choose correct assembly and choose replicable peaks
            # if want a specific replicate then change i["biological_replicates"] == 'rep_1' or 'rep_2'
            if len(i["biological_replicates"]) == 2 and i['assembly'] == assembly:
                bed_accession.append(i['accession'])

                
    # download all IDR bed files in accession
    for acc in bed_accession:
        download_url = f'https://www.encodeproject.org/files/{acc}/@@download/{acc}.bed.gz'
        bedgzfile = requests.get(download_url)
        if outpath:
            outfile = outpath + f'/{acc}.bed.gz'
        else:
            outfile = f'/{acc}.bed.gz'
        
        open(outfile, 'wb').write(bedgzfile.content)
    
    print(f'Download finished. All files are stored in {outpath} dir.')
            
     
def parseArgs():
    prs = argparse.ArgumentParser()
    prs.add_argument("version", type=str, help='Choose reference genome version (hg19 or GRCh38).')
    prs.add_argument("outdir",type=str, help="Path to directory to download file in")
    args = prs.parse_args()
    return args
    
if __name__=="__main__":
    args = parseArgs()
    getdata(args.version, args.outdir)

