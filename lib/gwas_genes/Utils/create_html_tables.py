import pandas as pd
import logging
import os

logging.basicConfig(format='%(created)s %(levelname)s: %(message)s', level=logging.INFO)

def create_datatable_html(csv_filepath, output_dir, rows_per_page=10):
    """
    Reads a CSV file, generates an HTML page with an interactive DataTable using DataTables.js,
    and saves it to the specified output directory.

    Args:
        csv_filepath (str): Path to the input CSV file.
        output_dir (str): Directory to save the generated HTML file.
        rows_per_page (int): Number of rows to display per page in the DataTable (default: 10).
    """
    try:
        df = pd.read_csv(csv_filepath)
        df.fillna("", inplace=True)
    except FileNotFoundError:
        print(f"Error: CSV file not found at {csv_filepath}")
        return
    except pd.errors.ParserError as e:
        print(f"Error: Could not parse CSV file at {csv_filepath}: {e}")
        return
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    file_name = os.path.basename(csv_filepath).split('.')[0]
    output_filepath = os.path.join(output_dir, f"{file_name}_datatable.html")

    # HTML template with placeholders for the table data
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Interactive DataTable</title>
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.8/css/jquery.dataTables.min.css">
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/buttons/2.4.2/css/buttons.dataTables.min.css">
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/select/1.7.0/css/select.dataTables.min.css">
        <style>
            body {{
                font-family: sans-serif;
            }}
            #datatable-container {{
                width: 95%;
                margin: 20px auto;
            }}
        </style>
    </head>
    <body>
        <div id="datatable-container">
            <table id="datatable" class="display compact" style="width:100%">
                <thead>
                    <tr>
                        {''.join(f'<th>{col}</th>' for col in df.columns)}
                    </tr>
                </thead>
                <tbody>
                    {''.join(f'<tr>{"".join(f"<td>{row[col]}</td>" for col in df.columns)}</tr>' for _, row in df.iterrows())}
                </tbody>
                <tfoot>
                    <tr>
                        {''.join(f'<th>{col}</th>' for col in df.columns)}
                    </tr>
                </tfoot>
            </table>
        </div>

        <script src="https://code.jquery.com/jquery-3.7.0.js"></script>
        <script src="https://cdn.datatables.net/1.13.8/js/jquery.dataTables.min.js"></script>
        <script src="https://cdn.datatables.net/buttons/2.4.2/js/dataTables.buttons.min.js"></script>
        <script src="https://cdn.datatables.net/select/1.7.0/js/dataTables.select.min.js"></script>
        <script>
            $(document).ready(function() {{
                $('#datatable').DataTable( {{
                    dom: 'Bfrtip',
                    buttons: [
                        'colvis', 'copy', 'excel', 'csv', 'pdf', 'print'
                    ],
                    select: true,
                    pageLength: {rows_per_page},
                    lengthMenu: [ 10, 25, 50, 100 ],
                    scrollX: true,
                    initComplete: function () {{
                        this.api()
                            .columns()
                            .every(function () {{
                                let column = this;
                                let title = column.header().textContent;
                                let input = document.createElement('input');
                                input.placeholder = title;
                                column.footer().replaceChildren(input);

                                input.addEventListener('keyup', () => {{
                                    if (column.search() !== this.value) {{
                                        column.search(input.value).draw();
                                    }}
                                }});
                            }});
                    }},
                }} );
            }} );
        </script>
    </body>
    </html>
    """

    # Write the HTML to the output file
    try:
        with open(output_filepath, "w") as f:
            f.write(html)
        print(f"DataTable HTML saved to {output_filepath}")
    except Exception as e:
        print(f"Error writing to file: {e}")


def create_index_page(csv_files, output_dir, output_filepath="index.html"):
    """
    Creates an index HTML page with tabs linking to DataTables for each CSV file.

    Args:
        csv_files (list): List of paths to CSV files.
        output_dir (str): Directory where the HTML files are stored.
        output_filepath (str): Path to save the generated index HTML file (default: "index.html").
    """
    index_html = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>CSV DataTables Index</title>
        <link rel="stylesheet" href="https://code.jquery.com/ui/1.12.1/themes/base/jquery-ui.css">
        <style>
            body {
                font-family: Arial, sans-serif;
                margin: 20px;
            }
            .tabs {
                margin-top: 20px;
            }
            .tabs ul {
                list-style: none;
                padding: 0;
            }
            .tabs li {
                display: inline-block;
                margin-right: 10px;
            }
            .tabs a {
                text-decoration: none;
                padding: 10px 20px;
                background-color: #f4f4f4;
                border-radius: 5px;
                color: #333;
                font-weight: bold;
            }
            .tabs a:hover {
                background-color: #ddd;
            }
        </style>
    </head>
    <body>
        <h1>CSV DataTables Index</h1>
        <div class="tabs">
            <ul>
    """
    
    # Generate tab links for each CSV file
    for csv_file in csv_files:
        file_name = os.path.basename(csv_file).split('.')[0]
        index_html += f'<li><a href="{output_dir}/{file_name}_datatable.html">{file_name}</a></li>\n'

    # Closing the tabs section and HTML
    index_html += """
            </ul>
        </div>
    </body>
    <script src="https://code.jquery.com/jquery-3.7.0.js"></script>
    <script src="https://code.jquery.com/ui/1.12.1/jquery-ui.js"></script>
    <script>
        $(document).ready(function() {
            $(".tabs").tabs();
        });
    </script>
    </html>
    """
    
    # Write the index HTML to the file
    try:
        with open(output_filepath, "w") as f:
            f.write(index_html)
        print(f"Index page saved to {output_filepath}")
    except Exception as e:
        print(f"Error writing index file: {e}")


if __name__ == "__main__":
    csv_files = ["file1.csv", "file2.csv", "file3.csv"]  # List your CSV files here
    csv_files = ["../../../../gwas_genes/test_local/workdir/tmp/gwas_genes_output/gwas_analysis_0_gene_centric_GWAS_d10000_p1e-05.csv","../../../../gwas_genes/test_local/workdir/tmp/gwas_genes_output/gwas_analysis_1_gene_centric_GWAS_d10000_p1e-05.csv"
    ]

    output_dir = "output_html"  # Specify your output directory here
    for csv_file in csv_files:
        create_datatable_html(csv_file, output_dir, rows_per_page=50)
    
    create_index_page(csv_files, output_dir, "output_html/index.html")

