import uuid
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport

class HTMLReportCreator:
    def __init__(self, callback_url):
        """
        Initializes the HTMLReportCreator class with the required clients.
        :param callback_url: The callback URL for KBase clients.
        """
        self.callback_url = callback_url
        self.dfu = DataFileUtil(callback_url)
        self.report = KBaseReport(callback_url)

    def create_html_report(self, output_dir, workspace_name, objects_created, all_results_path):
        """
        Creates an HTML report and uploads it to the KBase workspace.
        :param output_dir: Path to the directory containing the report files.
        :param workspace_name: Name of the workspace where the report will be stored.
        :return: Dictionary with the report name and reference.
        """
        report_name = 'gwas_genes_' + str(uuid.uuid4())
        file_links = []
        file_links.append({
                'path': all_results_path,
                'name': 'results_csv.zip'
        })


        # Upload the directory to Shock
        report_shock_id = self.dfu.file_to_shock({'file_path': output_dir, 'pack': 'zip'})['shock_id']

        # Create the HTML file metadata
        html_file = {
            'shock_id': report_shock_id,
            'name': 'index.html',
            'label': 'index.html',
            'description': 'HTML report for gwas_genes'
        }

        # Create the extended report
        report_info = self.report.create_extended_report({
            'direct_html_link_index': 0,
            'html_links': [html_file],
            'file_links': file_links,
            'report_object_name': report_name,
            'objects_created': objects_created,
            'workspace_name': workspace_name
        })

        return {
            'report_name': report_info['name'],
            'report_ref': report_info['ref']
        }

