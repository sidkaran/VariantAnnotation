#############################################################
#######Variant Annotator Example Prototype####################
#######Sid Kamalakaran 06/20/2019 #######################
#############################################################
#A basic script that takes in a VCF, pulls annotations from Ensembl Variant Effect Predictor and ExAC
#and outputs a tab delimited text file 
#############################################################

## We are using the PyVCF module which iterates over VCF and can grab required fields easily w/o rewriting a variant class/method
from __future__ import print_function
import json
import requests
import sys
import vcf
import argparse

# A generic request function to get batch queries with RESTful interfaces
def postFromREST(serverID,headers,data):
	r = requests.post(serverID, headers=headers, data=data)
	print(r)
	if not r.ok:
		r.raise_for_status()
		sys.exit()
	decoded = r.json()
	return(decoded)

# build an array with ExAC formated variants, query REST server; returns a Json dump
def get_exac(var_array):
	server="http://exac.hms.harvard.edu"
	ext="/rest/bulk/variant"
	headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
	datajson=json.dumps(var_array)
	decoded=postFromREST(serverID=server+ext,headers=headers,data=datajson)
	return(decoded)

# build an hash with Ensembl formated variants, query VEP REST server; returns a Json dump
def get_vep(var_array):
	test=var_array
	print("vdsfds")
	datajson = json.dumps({ "hgvs_notations" : test })
	#server = "https://rest.ensembl.org"
	headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
	server = "https://grch37.rest.ensembl.org"
	ext = "/vep/human/hgvs"
	decoded=postFromREST(serverID=server+ext,headers=headers,data=datajson)
	return(decoded)
### Main module to process the VCF and query/format json
def main():
	exacVars=[] # Initialize ExAC formated variant list
	ensemblVars=[] # Initialize Ensembl formated variant list
	maxLen=200 # determined max numer of variants in 1 query
	var_out={}
	vcf_reader = vcf.Reader(open(vcfFil, 'r'))
	for record in vcf_reader:
		for i in range(len(record.ALT)):# if multiple alt alleles are present, we expand the variant list to get individual effects
			#print(str(i))
			exacVars.append(str(record.CHROM) + "-" + str(record.POS)+"-"+record.REF+"-"+str(record.ALT[i]))
			ensemblVars.append(str(record.CHROM) + ":g." + str(record.POS)+record.REF+">"+str(record.ALT[i]))
			var_out[str(record.CHROM) + ":g." + str(record.POS)+record.REF+">"+str(record.ALT[i])]=[[],[],[]]
			var_out[str(record.CHROM) + ":g." + str(record.POS)+record.REF+">"+str(record.ALT[i])][0]=[str(record.CHROM),str(record.POS),str(record.REF),str(record.ALT[i]),str(record.INFO['DP']),str(record.genotype('vaf5')['DPR'][0]),str(record.genotype('vaf5')['DPR'][i+1]),str(record.genotype('vaf5')['DPR'][i+1]/(record.genotype('vaf5')['DPR'][0]+record.genotype('vaf5')['DPR'][i+1]))]


	#output a few examples for debugging
	print(exacVars[1:10])
	print(ensemblVars[1:10])	
	print(len(exacVars),len(ensemblVars),len(var_out.keys()))
	#exacVars=exacVars[1:70]
	#ensemblVars=ensemblVars[1:70]
	for i in range(0, len(ensemblVars),maxLen):
		ensemblL = ensemblVars[i:i + maxLen]
		#Grab chunks of 200 variants to requests
		decoded=get_vep(ensemblL)
		print("chunk"+str(i)+" retrieved")
		#decoded=get_vep(ensemblVars[1:2])
		for dr in decoded:
			ensemblAnnot=[dr["most_severe_consequence"],[]]
			#if dr.has_key("transcript_consequences"):
			if "transcript_consequences" in dr:
				for drc in dr["transcript_consequences"]:
					print(dr["input"],drc["consequence_terms"],drc["impact"],drc["gene_id"],drc["transcript_id"])
					if dr["most_severe_consequence"] in drc["consequence_terms"]:
						ensemblAnnot[1].append(drc["transcript_id"])
			else:
				ensemblAnnot[1].append("None")
			var_out[dr["input"]][1]=[ensemblAnnot[0],"|".join(ensemblAnnot[1])]


	print("Exac variants")

	for i in range(0, len(exacVars),maxLen):
		exacL = exacVars[i:i + maxLen]
		##Grab chunks of 200 variants to requests
		decoded=get_exac(exacL)
		print("chunk"+str(i)+" retrieved")
		for k in decoded:
			if "allele_freq" in decoded[k]["variant"]:
				#print(k,decoded[k]["variant"]["allele_freq"],decoded[k]["variant"]["allele_count"],decoded[k]["base_coverage"][0]["mean"])
				exacAnnot=[decoded[k]["variant"]["allele_freq"],decoded[k]["variant"]["allele_count"],decoded[k]["base_coverage"][0]["mean"]]
				varSplit=k.split("-")
				#print(varSplit)
				var_out[str(varSplit[0])+":g."+str(varSplit[1])+str(varSplit[2])+">"+str(varSplit[3])][2]=exacAnnot

	print(len(var_out.keys())) # tallyup total variants counts
	with open(outFil,'w') as out:
		out.write("VarID\tChr\tPos\tRef\tAlt\tReadDepthAtLocus\tRefReads\tAltReads\tVarAlleleFreq\tMostDamagingEffect\tTranscriptId(s)\tExAC AF\tExACAlleleCounts\tExACMeanCov\n")
		ctr=1
		for vo in var_out.keys():
			#print(vo,var_out[vo]); add filler elements when no info is retrieved from servers
			if var_out[vo][1]==[]:
				var_out[vo][1]=["None","None"]
			if var_out[vo][2]==[]:
				var_out[vo][2]=["NA","NA","NA"]
			out.write(vo+"\t"+"\t".join(str(i) for i in var_out[vo][0])+"\t"+"\t".join(str(i) for i in var_out[vo][1])+"\t"+"\t".join(str(i) for i in var_out[vo][2])+"\n")
			ctr+=1
		print(ctr,len(var_out.keys()))



if __name__== "__main__":
	parser = argparse.ArgumentParser(description='Annotate variants in VCF format')
	parser.add_argument('vcfFile', type=str, help='Input VCF file in 4.0 or above spec')
	parser.add_argument('tsvOut', type=str, help='Annotated TSV file')
	args = parser.parse_args()
	print(args.vcfFile)
	print(args.tsvOut)
	vcfFil=args.vcfFile
	outFil=args.tsvOut

	main()


