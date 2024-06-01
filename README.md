# DAESC-GPU
DAESC-GPU: A GPU-powered scalable software for single-cell allele-specific expression analysis

Here are the codes of DAESC-GPU. We use GPU to accelerate the excution of the statistical model DAESC. 

## Notebooks:
* DAESC_Small.ipynb : The notebook is used to fit true genetics data (400 genes, 30K+ cells). Details of running time can be found there. Generally speaking, the running time is 128s (T4), 71s (A100). Type I errors are 0.0675 (T4), 0.0725 (A100).
* DAESC_Large.ipynb : The notebook is used to fit true genetics data (4173 genes, 30K+ cells). Details of running time can be found there. Generally speaking, the running time is 570s (T4), 137s (A100).

## Usage: 
We heavily recommend one to use Google Colab to run the notebook above. Because we use the platform to develop the sofware. Additionally, anyone is available to free T4 GPU on Google Colab. You only need to change the address to import your data. Your data should be specific type. 

If one wants to repeat our work, you can send me an email to get the acsess to my processed data.  (tfcui23@uw.edu)
