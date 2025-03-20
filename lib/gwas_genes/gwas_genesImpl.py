# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import json
import csv
import pandas as pd

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil
from .Utils.get_gene_function import analyze_snps_and_genes
from .Utils.create_html_tables import create_datatable_html, create_index_page
from .Utils.html_report_creator import HTMLReportCreator
from .Utils.zip_files import copy_and_zip_csvs
#END_HEADER


class gwas_genes:
    '''
    Module Name:
    gwas_genes

    Module Description:
    A KBase module: gwas_genes
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        #self.shared_folder = "/kb/module/work"
        self.ws_url = config['workspace-url']
        self.ws_client = Workspace(self.ws_url)
        self.dfu = DataFileUtil(self.callback_url)
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_gwas_genes(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_gwas_genes
        logging.info("Starting run_gwas_genes with params: " + str(params))
        
        # Create output directory
        output_dir = os.path.join(self.shared_folder, 'gwas_genes_output')
        os.makedirs(output_dir, exist_ok=True)
        
        # Write parameters to file
        params_file = os.path.join(output_dir, 'params.json')
        with open(params_file, 'w') as f:
            f.write(json.dumps(params, indent=4))
        
        gwas_association_objects = params['gwas_association_objects']
        genome_ref = params['genome_ref']
        
        logging.info(f"Downloading Genome ref: {genome_ref}")
        # Get and write genome data
        genome_object_data = self.dfu.get_objects({'object_refs': [genome_ref]})['data'][0]['data']
        genome_file = os.path.join(output_dir, 'genome_data.json')
        with open(genome_file, 'w') as f:
            f.write(json.dumps(genome_object_data, indent=4))
        
        # Process each GWAS association object and save to separate files
        gwas_files = []
        for i, gwas_association_object in enumerate(gwas_association_objects):
            # Get the GWAS association object
            logging.info(f"Downloading GWAS association object: {gwas_association_object}")
            gwas_association_object_data = self.dfu.get_objects({'object_refs': [gwas_association_object]})['data'][0]['data']
            
            # Write each GWAS association object to a separate file
            gwas_file = os.path.join(output_dir, f'gwas_association_data_{i}.json')
            with open(gwas_file, 'w') as f:
                f.write(json.dumps(gwas_association_object_data, indent=4))
            
            gwas_files.append(gwas_file)

        logging.info(f"GWAS files: {gwas_files}")

        
        
        # Set analysis parameters
        distance_threshold = 10000  # 10kb
        pvalue_threshold = 1E-5    # p < 0.001
        all_csvs = list(); 
        for i, gwas_file in enumerate(gwas_files):
            # Convert GWAS file to SNP file format expected by the analysis function
            # First, extract just the data part we need
            with open(gwas_file, 'r') as f:
                gwas_data = json.load(f)

            
            # Run the analysis for this GWAS file
            try:
                logging.info(f"Running analysis for {gwas_file}")
                result = analyze_snps_and_genes(
                    gene_file=genome_file,
                    snp_file=gwas_file,
                    distance_threshold=distance_threshold,
                    pvalue_threshold=pvalue_threshold,
                    output_prefix=os.path.join(output_dir, f'gwas_analysis_{i}'),
                    save_output=True,
                    verbose=True,
                    save_gene_function= False  # Only save gene function for first analysis
                )
                
                all_csvs.append(result['snp_analysis_file'])
                all_csvs.append(result['gene_centric_file'])
                
                # Log the summary
                logging.info(f"Analysis {i+1} for {gwas_association_objects[i]} complete")
                logging.info(f"Summary: {result['summary']}")
                
            except Exception as e:
                logging.error(f"Error analyzing GWAS data {i}: {str(e)}")
        
        csv_zip_dir = os.path.join(self.shared_folder, "results_csv")
        output_dir1 = os.path.join(self.shared_folder, "results")

        zip_path = copy_and_zip_csvs(all_csvs, csv_zip_dir, result_csvs.zip)

        os.makedirs(output_dir1, exist_ok=True)
        for csv_file in all_csvs:
            create_datatable_html(csv_file, output_dir1, rows_per_page=10)
        create_index_page(csv_files, output_dir1, "index.html")
        output = {}
        report_creator = HTMLReportCreator(self.callback_url)
        objects_created = []
        output = report_creator.create_html_report(output_dir1, workspace, objects_created, zip_path)
        logging.info (output)


        
        #output = {
        #    'report_name': report_info['name'],
        #    'report_ref': report_info['ref'],
       # }

        #END run_gwas_genes

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_gwas_genes return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
