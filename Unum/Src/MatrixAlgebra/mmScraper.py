from urllib.request import urlopen, Request, urlretrieve
from bs4 import BeautifulSoup as soup
import os, sys
import gzip, shutil

"""
Extracts all files from mm to given directory. 
"""
def getAllMatrices(directory):
    
    if not os.path.isdir(directory): os.mkdir(directory)

    base_url        = 'https://math.nist.gov/MatrixMarket/matrices.html'
    base_matrix_url = 'https://math.nist.gov'

    req = Request(base_url)
    with urlopen(req) as page:
        matrices = soup(page.read(), 'html.parser')

    table          = matrices.find('table')
    matrix_columns = table.findAll('td'   )

    matrix_links = [] #store all the links here
    for col in matrix_columns:
        matrix_col = col.findAll('a')
        for matrix in matrix_col:
            matrix_links.append(base_matrix_url + matrix['href'])

    download_links = []
    for matrix in matrix_links:
        req = Request(matrix)
        with urlopen(req) as matrix_page:
            matrix_soup = soup(matrix_page.read(), 'html.parser')
            li = matrix_soup.find('li')
            download_links.append(li.findAll('a')[1]['href'])

    for dl in download_links:
        matrix_name  = dl.split('/')[-1].replace('.mtx.gz', '')
        compressed   = directory + '/' + matrix_name + '.mtx.gz'
        mtx_filename = directory + '/' + matrix_name + '.mtx'

        urlretrieve(dl, compressed)

        with gzip.open(compressed, 'rb') as f_in:
            with open(mtx_filename, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        os.remove(compressed)
        
if __name__ == "__main__":
    getAllMatrices(sys.argv[1])