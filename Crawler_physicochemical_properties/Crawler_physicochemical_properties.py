import pandas as pd
import time
import requests

# Fetch physical and chemical properties Property_Name
# Property_Name
Property_Name = ['Molecular Weight',
                'XLogP3',
                'Hydrogen Bond Donor Count',
                'Hydrogen Bond Acceptor Count',
                'Rotatable Bond Count',
                'Exact Mass',
                'Monoisotopic Mass',
                'Topological Polar Surface Area',
                'Heavy Atom Count',
                'Formal Charge',
                'Complexity',
                'Isotope Atom Count',
                'Defined Atom Stereocenter Count',
                'Undefined Atom Stereocenter Count',
                'Defined Bond Stereocenter Count',
                'Undefined Bond Stereocenter Count',
                'Covalently-Bonded Unit Count',
                'Compound Is Canonicalized',
                'Physical Description', 
                'Color/Form', 
                'Odor', 
                'Boiling Point',
                'Melting Point',
                'Flash Point', 
                'Solubility',
                'Density']

def get_Pubchem(url):
    '''
    Retrieve compound data from PubChem.
    Input:
        url: API address obtained from PubChem
            For example: https://pubchem.ncbi.nlm.nih.gov/compound/2244
            The address changes depending on the compound's CID (PubChem Compound ID).
    Output:
        JSON data returned by PubChem
    '''

    # The API address to be used
    # url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/2244/JSON/?heading=Names+and+Identifiers'

    try:
        # Send a request to fetch data
        response = requests.get(url=url, timeout=10)
        return response
    except requests.exceptions.Timeout:
        print("Request timed out. Please try again later.")
        return -1
    except requests.exceptions.ConnectionError as e:
        print(f"Connection error: {e}")
        return -1
    except requests.exceptions.RequestException as e:
        print(f"An exception occurred: {e}")
        return -1

# Parse name from response
def get_name_value(response, name):
    '''
    Input:
        response: The response object
        name: Attributes such as Canonical SMILES, Molecular Formula, CAS, Molecular Weight, etc.
    Output:
        Corresponding values of attributes like SMILES, MF, CAS, Molecular Weight, etc.
    '''
    # Use the split method for parsing. This minimizes errors.
    # Alternatively, convert the response to JSON and parse it, but this might raise errors for missing attributes.

    # If the specified name is not found, split will raise an error, and None will be returned in such cases.
    try:
        data = response.text
        data = data.split(f'"TOCHeading": "{name}",')[1].split('"Information": [')[1].split('},')[0]
        data = data.replace('\n ', '').replace(' ', '')
        try:
            # Parse Chemical and Physical Properties
            Value = data.split('String":"')[1].split('"}')[0]
        except:
            # Parse Chemical and Physical Properties (alternative method)
            Value = data.split('"Number":')[1].split('}')[0].replace('[', '').replace(']', '')
        return Value
    except:
        return None

# Parse multiple values for a given name
def get_name_value2(response, name):
    '''
    Input:
        response: The response object
        name: Attributes like Canonical SMILES, Molecular Formula, CAS, Molecular Weight, and others
    Output:
        A list of values corresponding to the attribute.
    '''
    # Parse the response using split, minimizing potential errors.
    # If the specified name is missing, None is returned.

    Value = []
    try:
        data = response.text
        data = data.split(f'"TOCHeading": "{name}",')[1].split('"Information": [')[1].split('},')
        data = [k.replace('\n ', '').replace(' ', '') for k in data]

        for k in data:
            if 'String":"' in k:
                v = k.split('String":')[1].split('}')[0]
                if '"Markup":' in v:
                    v = v.split(',"Markup"')[0]
                v = v.replace('\\u000', '')
                v = v[1:-1]
            elif '"Number":' in k:
                v = k.split('"Number":')[1].split('}')[0].replace('[', '').replace(']', '')
                v = v[1:-1]
            else:
                pass
            Value.append(v)

        Value = list(set(Value))
        return Value
    except:
        return []

# Read the CSV file
csv_file = 'your_csv_file.csv'
print(csv_file)
df = pd.read_csv(csv_file)

# Get the values from the "CID" column and convert them to a list
cid_list = df['CID'].tolist()

# Fetch data in batches for names, IUPAC Name, InChI, InChI Key, SMILES, Molecular Formula, CAS, and other properties based on CIDs
  
df = pd.DataFrame()

count = 0
for cid in cid_list[:]:
    if cid != 0:
        # API address for Names and Identifiers
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/'

        response = get_Pubchem(url)
        if response != -1:
            # Parse name and CID
            try:
                cid_paras = response.json()['Record']['RecordNumber']
                compound_name = response.json()['Record']['RecordTitle']
      
                d = {}
                d['cid'] = cid
                d['compound_name'] = compound_name
                d['cid_paras'] = cid_paras

                # Parse property names
                for name in Property_Name:
                    Value = get_name_value(response=response, name=name)
                    d[name.replace(' ', '_').replace('-', '_')] = Value

                df_tmp = pd.DataFrame([d])
                df = pd.concat([df, df_tmp])

                print(f'cid: {cid}, status_code: {response.status_code}, df: {df.shape[0]}')
            except KeyError as e:
                print(f"KeyError occurred: {e}. Skipping.")
            except requests.exceptions.RequestException as e:
                print(f"An exception occurred: {e}")
    else:
        d = {}

    # Pause for 0.5 seconds
    time.sleep(0.5)
    count += 1
    if count % 10 == 0:
        time.sleep(2)

df = df.reset_index(drop=True)

# Write the DataFrame to a CSV file
output_file = 'your_output_file.csv'
df.to_csv(output_file, index=False)
print("DataFrame successfully written to CSV file.")
