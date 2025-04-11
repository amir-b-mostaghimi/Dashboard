import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
from stmol import showmol

# Load the dataframes
all_pairs_df = pd.read_csv('/Users/amir/Dropbox/Set_up/REE_pairs.csv')

df = pd.read_csv('/Users/amir/Dropbox/Set_up/REE_D.csv')

def find_reference(df, target_smiles):
    target_idx = df[df['SMILES'] == target_smiles].index[0]
    nan_rows_above = df.loc[:target_idx].loc[df['REE'].isna()]
    if not nan_rows_above.empty:
        return nan_rows_above.iloc[-1]
    return None



def makeblock(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mblock = Chem.MolToMolBlock(mol)
    return mblock

def render_mol(xyz):
    xyzview = py3Dmol.view()
    xyzview.addModel(xyz, 'mol')
    xyzview.setStyle({'stick':{}})
    xyzview.setBackgroundColor('white')
    xyzview.zoomTo()
    showmol(xyzview, height=500, width=500)

def create_dashboard(all_pairs_df, df):
    st.title("Lanthanide Separation Factors Dashboard")
    
    # Add toggle for 3D visualization
    show_3d = st.sidebar.checkbox("Show 3D Molecular Structures", value=False)
    
    # Get unique first elements from the Pair column
    first_elements = sorted(list(set([pair.split('/')[0] for pair in all_pairs_df['Pair'].unique()])))
    
    # First dropdown for selecting first element
    first_element = st.selectbox(
        'Select first lanthanide:',
        first_elements
    )
    
    # Get available second elements for the selected first element
    available_seconds = sorted(list(set([
        pair.split('/')[1] 
        for pair in all_pairs_df['Pair'].unique() 
        if pair.split('/')[0] == first_element
    ])))
    
    # Second dropdown for selecting second element
    second_element = st.selectbox(
        'Select second lanthanide:',
        available_seconds
    )
    
    # Create the pair string
    selected_pair = f"{first_element}/{second_element}"
    
    # Filter dataframe for selected pair
    pair_data = all_pairs_df[all_pairs_df['Pair'] == selected_pair].nlargest(3, 'Separation_Factor')
    
    if not pair_data.empty:
        st.header(f"Top 3 Separation Conditions for {selected_pair}")
        
        for i, row in pair_data.iterrows():
            st.subheader(f"Rank {i+1} (Separation Factor: {row['Separation_Factor']:.2f})")
            
            # Create two columns for better layout
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("**Reaction Conditions:**")
                st.write(f"- Solvent: {row['Solvent']}")
                st.write(f"- Acid Type: {row['Acid Type']}")
                st.write(f"- Acid Concentration: {row['Acid Conc (M)']} M")
                st.write(f"- Temperature: {row['Temperature (C)']}°C")
                st.write(f"- Ligand Concentration: {row['Ligand Conc (M)']} M")
            
            with col2:
                st.markdown("**Distribution Coefficients:**")
                st.write(f"- {first_element} D value: {row['Current_REE_D']:.2e}")
                st.write(f"- {second_element} D value: {row['Next_REE_D']:.2e}")
            
            st.markdown("**SMILES Structure:**")
            st.code(row['SMILES'], language='text')
            
            # Add reference information
            try:
                reference = find_reference(df, row['SMILES'])['SMILES']
                if reference is not None:
                    st.markdown("**Reference:**")
                    st.write(reference)  # Simply write the string directly
            except Exception as e:
                st.warning(f"Could not find reference for this molecule. Error: {str(e)}")
            # Add 3D visualization if enabled
            if show_3d:
                try:
                    st.markdown("**3D Structure:**")
                    xyz = makeblock(row['SMILES'])
                    render_mol(xyz)
                except Exception as e:
                    st.warning(f"Could not generate 3D structure for this molecule. Error: {str(e)}")
            
            st.markdown("---")
    else:
        st.warning("No data available for this pair")
    
    # Add some basic statistics
    st.sidebar.header("Dataset Statistics")
    st.sidebar.write(f"Total number of pairs: {len(all_pairs_df['Pair'].unique())}")
    st.sidebar.write(f"Number of unique SMILES: {len(all_pairs_df['SMILES'].unique())}")
    
    # Add range of separation factors
    st.sidebar.write(f"Separation Factor Range:")
    st.sidebar.write(f"Min: {all_pairs_df['Separation_Factor'].min():.2f}")
    st.sidebar.write(f"Max: {all_pairs_df['Separation_Factor'].max():.2f}")

if __name__ == "__main__":
    create_dashboard(all_pairs_df, df)
