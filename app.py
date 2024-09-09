import io
import matplotlib.pyplot as plt
from flask import Flask, render_template, request, jsonify
import requests
import json
from matplotlib.patches import Rectangle
from io import BytesIO
import base64
import matplotlib
import random
matplotlib.use('Agg')  # Use the non-GUI Agg backend
app = Flask(__name__)

# Fetch gene information from Ensembl
selected_variant = None

def get_gene_info(gene_name):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/homo_sapiens/{gene_name}?expand=1"
    headers = {"Content-Type": "application/json"}
    response = requests.get(server + ext, headers=headers)
    if not response.ok:
        return None
    return response.json()

# Fetch gene sequence from Ensembl


def get_gene_sequence(gene_id):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{gene_id}?"
    headers = {"Content-Type": "application/json"}
    response = requests.get(server + ext, headers=headers)
    if not response.ok:
        return None
    return response.json().get('seq')


def get_variants_in_gene(gene_id):
    server = "https://rest.ensembl.org"
    ext = f"/overlap/id/{gene_id}?feature=variation"
    headers = {"Content-Type": "application/json"}
    response = requests.get(server + ext, headers=headers)
    if not response.ok:
        # Return None for both variants and the selected variant if the API call fails
        return None, None
    variants = response.json()
    if not variants:
        return None, None
    selected_variant = random.choice(variants)  # Randomly select one variant
    return variants, selected_variant

# Plot exon locations and return as a base64 image

def plot_exon_locations(gene_info, gene_sequence):
    if not gene_info or not gene_sequence:
        return None

    plt.figure(figsize=(10, 2))
    ax = plt.gca()
    ax.add_patch(Rectangle((0, 0.5), len(gene_sequence), 0.5,
                 edgecolor='black', facecolor='lightgray'))
    for transcript in gene_info.get('Transcript', []):
        for exon in transcript.get('Exon', []):
            exon_start = exon['start'] - gene_info['start']
            exon_length = exon['end'] - exon['start'] + 1
            ax.add_patch(Rectangle((exon_start, 0.5), exon_length,
                         0.5, edgecolor='blue', facecolor='pink'))

    ax.set_xlim(0, len(gene_sequence))
    ax.set_ylim(0, 1)
    ax.set_xlabel('Nucleotide Position')
    plt.title(f"Exon Locations in Gene: {gene_info.get('display_name')}")

    # Convert plot to PNG image and encode it
    img = io.BytesIO()
    plt.savefig(img, format='png', bbox_inches='tight')
    plt.close()
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode('utf-8')
    return plot_url


def plot_variants(gene_info, variants, selected_variant):
    if not gene_info or not variants:
        return None

    gene_length = gene_info['end'] - gene_info['start'] + 1
    fig, ax = plt.subplots(figsize=(10, 2))
    ax.add_patch(Rectangle((0, 0.5), gene_length, 0.5,
                 edgecolor='black', facecolor='lightgray'))

    for variant in variants:
        variant_position = variant['start'] - gene_info['start']
        ax.plot([variant_position, variant_position],
                [0.4, 0.6], color='red', linewidth=2)

    ax.set_xlim(0, gene_length)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Nucleotide Position')
    ax.set_title(f"Variant {selected_variant.get('assembly_name')} in Gene: {
                 gene_info.get('display_name')}")

    # Convert plot to PNG image
    img = BytesIO()
    plt.savefig(img, format='png', bbox_inches='tight')
    plt.close(fig)
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode('utf-8')
    return plot_url


@app.route('/', methods=['GET', 'POST'])
def index():
    gene_info = None
    plot_url = None
    variant_plot_url = None 
    selected_variant = None  # Initialize the selected variant
    error = None

    if request.method == 'POST':
        gene_name = request.form['gene_name']
        gene_info = get_gene_info(gene_name)
        if gene_info:
            gene_sequence = get_gene_sequence(gene_info['id'])
            if gene_sequence:
                plot_url = plot_exon_locations(gene_info, gene_sequence)
                variants, selected_variant = get_variants_in_gene(
                    gene_info['id'])
                if variants:
                    variant_plot_url = plot_variants(
                        gene_info, variants, selected_variant)
                else:
                    error = "No variants found or API error."
            else:
                error = "Could not retrieve gene sequence."
        else:
            error = "Gene not found or API error."

    return render_template('index.html', gene_info=gene_info, plot_url=plot_url, variant_plot_url=variant_plot_url, selected_variant=selected_variant, error=error)


if __name__ == '__main__':
    app.run(debug=True)
