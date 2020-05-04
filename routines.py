# Authors: Aris Marcolongo & Tobias Binninger
# Testing in progress ( 03 May 2020 ) : Steven D. Lacey & Federico Zipoli
# Affiliation: IBM

#
# Routines to compute electrochemical stability using the phase stability method as described in :  
# "Comparison of computational methods for the electrochemical stability window of solid-state electrolyte materials"
# Journal of Materials Chemistry A, 3 (2020) . This method is an alternative to the standard convex-hull analyisis
# and focuses on the instability reactions driving the decomposition of a selected SSE, after contact with an electrode
# material.
# This version uses data downloaded from the Materials Project
#
# Contains:
#
# - stability_window : main driver
# - download_and_filter_types, download_and_filter_formula : functions used to download the needed data from MP
# - filter_compounds_formula : it filter the compunds downloaded according to requested critera and evaluating minimum energy compounds
# - expand_stoichiometry : checks if a reaction is stoichiometrically fisible
# - compute_potential : evaluates potential at which the reaction proceeds, with an error estimate
# - search_reagent : queries the MP project for the SSE under consideration

import numpy as np
import itertools
import time
from itertools import chain, combinations
from pymatgen.core.composition import Composition


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
    

def download_and_filter_types(types, mpr, icsd, verbosity=True, remove_non_comp_ox = -1  ):
     
        products=[] 
        
        if (verbosity):
            print('download_and_filter_types')
            
        ntypes=len(types)
 
        tic=time.time()
        for item in powerset(types):
            if (len(item)>0):
                if (verbosity):
                    print(item)
                data=mpr.query(criteria={"elements":{"$all": list(item)},"nelements": len(item)},
                     properties=['formula','final_energy','icsd_ids','spacegroup','energy_per_atom',
                                 'nsites','volume','material_id','pretty_formula'])
                products.append(data)
                if (verbosity):
                    print('Elapsed time for downloading: ',time.time()-tic)
                    
        print('Elapsed time for downloading: ',time.time()-tic)
        
        if (verbosity):
            print('Made queries: ',len(products))
        
        formulas_downloaded=[]
        products_downloaded={}
        
        for iset in products:
            for iproduct in iset:
                if (iproduct['formula'] not in formulas_downloaded):
                    formulas_downloaded.append(iproduct['formula'])
                    products_downloaded[len(formulas_downloaded)-1]=[]
                    products_downloaded[len(formulas_downloaded)-1].append(iproduct)
                else:
                    n=formulas_downloaded.index(iproduct['formula'])
                    products_downloaded[n].append(iproduct)
        
        products_cleaned=[]
        for iprod in range(len(products_downloaded)):
            
            filtered=filter_compounds_formula(products_downloaded[iprod], verbosity=verbosity, 
                                              icsd=icsd, remove_non_comp_ox = remove_non_comp_ox )
            
            #filtered can be None if there are no structures in products_downloaded[iprod] with an ICSD entry (when icsd = True)
            if (filtered is not None):
                products_cleaned.append(filtered)
            else:
                print('Not in ICSD: ', products_downloaded[iprod][0]['formula'])
                
        if (verbosity):
            print(' End download and filter types---- \n')
            
        return products_cleaned
    
    
def download_and_filter_formula(formulas, mpr, icsd, verbosity=True, remove_non_comp_ox = -1 ):
     
        if (verbosity):
            print('download_and_filter_formula')
            
        products=[]
        for formula in formulas:

                print(formula)
                   
                tic=time.time()
                while (True):
            
                   try:
                    data=mpr.query(criteria={"formula": formula},
                               properties=['formula','final_energy','icsd_ids','spacegroup',
                                           'energy_per_atom','nsites','volume','material_id','pretty_formula'])
                    break
                
                   except Exception :
                    print('MPRestError --retry')
                    pass
            
                print('Elapsed time for downloading: ',time.time()-tic)  

                if (len(data)==0):
                    print('Nothing found')
                    return None
                
                data_chosen=filter_compounds_formula(data, verbosity=verbosity, icsd=icsd, remove_non_comp_ox = remove_non_comp_ox )
                
                 #data_chosen can be None if there are no structures in products_downloaded[iprod] with an ICSD entry
                if (data_chosen is not None):
                    products.append(data_chosen)
                else:
                    print('Not in ICSD: ', data[0]['formula'])
                
        if (verbosity):
            print(' End download and filter formula---- \n')
            
        return products  
    
    
def filter_compounds_formula(data, icsd ,verbosity=True, remove_non_comp_ox = -1 ):

    #check if at least one compound belongs to ICSD
    t_icsd=False
    for idata in range(len(data)):
        if len(data[idata]['icsd_ids'])>0 :
            t_icsd=True 
            break
            
    if (icsd):
        if (not t_icsd):
            return None
            
    if (verbosity):
            print('filter_compounds_formula')
            print(data[0]['formula'])
# Take only simulations with a certain number of sites "nsites_cutoff"

    nsite_max=data[0]['nsites']
    isite_max=0
    nsites=[] 
    for idata in range(len(data)):
        nsites.append(data[idata]['nsites'])
        if (data[idata]['nsites'] >= nsite_max):
            nsite_max=data[idata]['nsites']
            isite_max=idata
    if (verbosity):
        print(nsite_max,isite_max)

    nsites_cutoff=np.mean(np.array(nsites))

# Select simulations with the minimum energy     

    minimo=data[isite_max]['energy_per_atom']
    imin=isite_max
    cont=0
    for idata in range(len(data)):

        if (verbosity):
            print(idata,data[cont]['final_energy']/data[cont]['nsites'],data[cont]['energy_per_atom'],data[cont]['nsites'])

        if ( (data[cont]['energy_per_atom']<minimo) and (data[cont]['nsites']>nsites_cutoff) ):
            imin=cont
        cont+=1

    if (verbosity):
        print(imin)

    if (verbosity):
            print('End filter compounds formula---- \n')
            
    if ( remove_non_comp_ox > -1 ):
        
        # for pure elements we have no issue
        if ( len(data[imin]['formula'].keys()) > 1 )  :
    
            if (t_icsd):
                aos=True
            else:
                aos=False
                
            comp = Composition(data[imin]['pretty_formula']) 

            # if pymatgen does not find good oxidation states we have an issue
            if (len(comp.oxi_state_guesses(all_oxi_states=aos))==0 ):

                    if ( remove_non_comp_ox == 0 ):

                        print(data[imin]['pretty_formula'],' non compatible. Asked to remove')
                        return None

                    else:

                        found=False 

                        for nhyd in range(100):
                            app='H'+str(nhyd)
                            comp = Composition(app+data[imin]['pretty_formula'])
                                         
                            #when adding H we try to be more sure of the oxi state and we have therefore all_oxi_states = False
                            if (len(comp.oxi_state_guesses(all_oxi_states=False))>0 ):
                                print('aos ',aos)
                                found=True
                                if ('H' in data[imin]['formula'].keys()):
                                    data[imin]['formula']['H']=data[imin]['formula']['H']+nhyd
                                else:
                                    data[imin]['formula']['H']=nhyd
                                print(data[imin]['pretty_formula'],' non compatible. Asked to add ', str(nhyd), 'hydrogen atoms')
                                print('OX states: ', comp.oxi_state_guesses(all_oxi_states=False))
                                break

                        if not found :
                            print(data[imin]['pretty_formula'],' non compatible but did not find hydrogens to add.')
                            return None
                
    return data[imin]

def expand_stoichiometry(reagent, products, det_min=0, verbosity=False):
    
    """
  
    Parameters:
    
    The dictionary of the reagent and of the products containing their field formula.

    Returns:
    
    The stoichiometric coefficients balancing the reaction, if possible in the given order. If 
    not possible it returns None. Usually Lithium 
    should be the last product so that v_stoi_out[-1] > or < 0 indicates an oxidation or reduction reaction
    
    
    """
    if (verbosity):
            print('Expand Stoichiometry')
            
    reagent_elements = list(reagent[0]['formula'].keys())
    
    num_elements = len(reagent_elements)
    num_products = len(products)
 
    if (num_products > num_elements):
        print('The number of products should be less or equal to the number of elements')
        raise InputError 
        
    if (verbosity): 
        print('Reagent elements: ' ,reagent_elements)
    
    v_stoi_in = [reagent[0]['formula'][element] for element in reagent_elements]
    v_stoi_out = np.zeros(num_products)
    
    if (verbosity):
        print('Reagent: ')
        print(v_stoi_in)
    
    h = np.zeros((num_elements, num_products))
    for i in range(num_elements):
        for j in range(num_products):
            try:
                h[i,j] = products[j]['formula'][reagent_elements[i]]
            except KeyError :
                pass
            
    if (verbosity):
        print('The matrix of product coefficients ( columnwise ):')
        print(h)
    
    # We see if we can balance the first n_products elements by inverting the upper square matrix of 
    # dimension n_products x n_products .
    # Than we need to check if the solution applies to all elements in the case n_elements > n_products
    
    if  np.abs(np.linalg.det(h[:num_products,:num_products])) > det_min:
        
            v_stoi_out = np.matmul(np.linalg.inv(h[:num_products,:num_products]), v_stoi_in[:num_products])
            
            check=np.zeros(num_elements)
            for iprod in range(num_products):
                check=check+v_stoi_out[iprod]*h[:,iprod]
                
            if (np.abs(np.sum(check-v_stoi_in)))>0.0001 :
                v_stoi_out=None
            else:
                if (verbosity):
                    print('System solved check: ',np.abs(np.sum(check-v_stoi_in)))
            
    else:
            v_stoi_out = None
    
    if (verbosity):
        print(v_stoi_out)
        
        if (v_stoi_out is None):
            print('WARNING: the linear system was not solvable')
            
    if (verbosity):
            print('End expand stoichiometry---- \n')

    return v_stoi_out

def compute_potential(reagent, products, det_min=0, perc_err_energy=0.01, verbosity=True):
    
    """
  
    Parameters:
    
    - reagent, products: The dictionary of the reagent and of the products containing their field formula and the energies.
    - det_min : a threshold for the determinant for the matrix of stoichiometric coefficients to be invertible
    - perc_err_energy : a default percentage error of perc_err_energy on the total energy of each material is supposed

    Returns:
    
    A tuple (reaction_type ( vitual values in [ox,red,decomposition], pot). 
    If the reaction is (decomposition) it returns type,None,None, None
    
    If there are problems with the stoichiometry (non-positive stoichimoetries for the products or the reaction cannot proceed) 
    it returns None, None , None , None

    """
    if (verbosity):
            print('Compute_potential')
            
    #check that last product is lithium
    nat_formula=np.sum([ products[-1]['formula'][key] for key in products[-1]['formula'].keys() ])
    if  ('Li' not in products[-1]['formula'].keys()):
        raise InputError
    else:
        if (products[-1]['formula']['Li']<nat_formula):
            raise InputError
            
    nproducts=len(products)
    
    # Find coefficients of the reaction
    coeffs = expand_stoichiometry(reagent, products, det_min=0, verbosity=verbosity)
   

    # check coefficients
    if (coeffs is None):
        if (verbosity):
            print('WARNING: stoichiometry non compatible')
        return None, None, None, None
    if any(coeffs[i] < -0.0001 for i in range(nproducts-1)):
        if (verbosity):
            print('WARNING: stoichiometry cannot proceed')
        return None, None, None, None       
    
    nprod=len(products)-1
    
    # Find reaction type
    epsilon=1e-6
    if coeffs[-1] > epsilon:
        reaction_type='Oxidation'
    elif coeffs[-1] < -epsilon:
        reaction_type='Reduction'
    else:
        reaction_type= 'Phase Separation'

    if (reaction_type=='Phase Separation'):
        if (verbosity):
            print('WARNING: Phase Separation reaction. No potential computed. The reaction would proceed \
                without electrode. Neglected.')
        
        return reaction_type, None, None, None
        
        
    # Compute normalized energies, i.e. per formula unit
    normalized_prod_energies=[]
    for product in products:
        
        nat_formula=np.sum([ product['formula'][key] for key in product['formula'].keys() ])
        energy=product['final_energy']/( product['nsites']/nat_formula )
        
        normalized_prod_energies.append(energy)
        print(nat_formula,energy)
        
    nat_formula=np.sum([ reagent[0]['formula'][key] for key in reagent[0]['formula'].keys() ])
    normalized_reag_energy=reagent[0]['final_energy']/( reagent[0]['nsites']/nat_formula )
    
    
    # Compute potential 
    # Phi = - ( ( Sum_Pr ( coeff(Pr) * E_Pr ) - E_reag ) / x - E_lit )
    
    pot=0
    x = - coeffs[-1]
    for iprod in range(nproducts-1):
        pot = pot + coeffs[iprod] * normalized_prod_energies[iprod]
    pot = pot - normalized_reag_energy
    pot = pot/x
    pot = pot - normalized_prod_energies[-1]
    pot = - pot
    
    # compute estimated arror potential :
    # Err_Phi = perc_err_energy * np.sqrt( ( Sum_Pr coeff²(Pr) E_pr² + E_reag² )/x² + E_lit² ) 
    
    err_pot=0
    for iprod in range(nproducts-1):
        err_pot = err_pot + ( coeffs[iprod]**2 ) * (normalized_prod_energies[iprod]**2)
    err_pot = err_pot + normalized_reag_energy**2   
    err_pot = err_pot / (x**2)
    err_pot = err_pot + ( normalized_prod_energies[-1] )**2
    err_pot = perc_err_energy * np.sqrt(err_pot)
    
    if (verbosity):
        print(pot)
 
    if (verbosity):
            print('End compute potential---- \n')
            
    return  reaction_type, pot, err_pot, coeffs

def search_reagent(reagent, mpr, verbosity=True ):
    
    if (verbosity):
            print(' search_reagent--- \n')
            
    nreag_elements=len(reagent)
    print(reagent.keys())
    
    for permutation in itertools.permutations(reagent.keys()):
        
        tic=time.time()
        search={}
        for ireag in range(nreag_elements):
            search[permutation[ireag]]=reagent[permutation[ireag]]
        print('search: ', search)
        while (True):
            
            try:
                data=mpr.query(criteria={"formula": search},
                   properties=['formula','final_energy','icsd_ids','spacegroup','energy_per_atom',
                               'nsites','volume'])
                break
                
            except Exception :
                
                print('MPRestError --retry')
                pass
            
        print('Elapsed time for downloading: ',time.time()-tic)       
        if (len(data)>0):
            print('compund_found')
            break
            
        if (verbosity):
            print(data)
            
    # In this case for the reagent we never ask to be part of ICSD        
    data_chosen=filter_compounds_formula(data, verbosity=verbosity, icsd=False, remove_non_comp_ox = -1 )        
    
    if (verbosity):
            print(' End search reagent---- \n')
                        
    return data

def stability_window(reagent, mpr, verbosity=True, icsd = False, remove_non_comp_ox = -1 ):
    
    """ Main driver to evaluate the electro-chemical stability window from MP data
  
    Parameters:
    
    - reagent: the dictionary of the reagent, i.e. the material for which the stability window needs to be evaluated.
    - mpr: an MPRester object used for the Materials Project Rest API, i.e. to download the data from the internet.
    - icsd : if True, only compounds appearing both on MP and on ICSD are considered
    - remove_non_comp_ox : integer variable assuming values in [-1,0,1]. If remove_non_comp_ox == -1, all compounds are kept,
                            if remove_non_comp_ox == 0 charged materials are remove,  if remove_non_comp_ox == 1 attached 
                            hydrogens are added to negatively charged materials. remove_non_comp_ox == 1 may be used only
                            for SSE containing hydrogen.

    Returns:
    
    - oxidations: a dictionary oxidations with all possible oxidation reactions. The keys  "pot", "error", "products" "coeffs" refer to the       potential, the estimated error, the products of the reaction
      and the stoichiometric coefficients.
    
    - reductions: the same for the possible reductions.

    """
        
    if (remove_non_comp_ox == -1):
        print('remove_non_comp_ox to -1 : we keep all compounds')
        
    elif (remove_non_comp_ox == 0):
        print('remove_non_comp_ox to 0 : we remove non stoichio compounds')
        
    elif (remove_non_comp_ox == 1):
        print('WARNING remove_non_comp_ox to +1 : we add hydrogens to non stoichio compounds. Non tested. May work if the SSE contains hydrogen but not in other cases at the moment. ')
        
    else : 
        print('Non admissible value for remove_non_comp_ox')
        raise Exception
        
    if (verbosity):
            print(' stability_window')
            
    reagent = search_reagent(reagent, mpr, verbosity=verbosity )
    nreag_elements=len(list(reagent[0]['formula'].keys()))
    print('nreag_elements', nreag_elements)
    
    if  ( ('H' not in list(reagent[0]['formula'].keys()) ) and (remove_non_comp_ox == 1) ):
        print('H not in reagent and remove_non_comp_ox == 1. STOP')
        raise Exception
    
    products = download_and_filter_types(list(reagent[0]['formula'].keys()), mpr, verbosity=verbosity, icsd=icsd,
                                         remove_non_comp_ox = remove_non_comp_ox )
    
    print('List of final possible products: ')
    for product in products:
        print(product['formula'])
        
        
    lithium=None
    for product in products:
        if (list(product['formula'].keys())==['Li']):
            print('Found Lithium in the products')
            lithium=product
            products.remove(product)
    if (lithium is None):
        print('ERROR LITHIUM NOT FOUND')
        raise Exception
    
    found=False
    for product in products:
        if (product['formula']==reagent[0]['formula']):
            print('Reagent found in products')
            found=True
            products.remove(product)      
    # if icsd the reagent may be not found because it is not on ICSD        
    if  ( (not found) and (not icsd) ):
        print('ERROR REAGENT NOT FOUND')
        raise Exception
        
    candidates=len(products)
    print('Candidate_decomposition_products: ', candidates)
    
    oxidations={}
    oxidations['pot']=[]
    oxidations['products']=[]
    oxidations['error']=[]
    oxidations['coeffs']=[]
    
    reductions={}
    reductions['pot']=[]
    reductions['products']=[]
    reductions['error']=[]
    reductions['coeffs']=[]
    
    
    for k in range(nreag_elements-1):
    
        for item in itertools.combinations(range(candidates),k+1):
 

                print('Hi reaction considered: ', item, k+1)
                products_=[]
                for jitem in item:
                    products_.append(products[jitem])
                products_.append(lithium)
 
                
                reaction_type,pot,err_pot, coeffs =compute_potential(reagent, products_, 
                                                                     det_min=0, verbosity=False)

                if ( (reaction_type is not None) and (pot is not None) ):

                    print('Reaction_found')
                    for product in products_:
                        print(product['formula'])

                    if (reaction_type=='Oxidation'):
                        oxidations['pot'].append(pot)
                        oxidations['products'].append(products_)
                        oxidations['error'].append(err_pot)
                        oxidations['coeffs'].append(coeffs)
                    

                    if (reaction_type=='Reduction'):
                        reductions['pot'].append(pot)
                        reductions['products'].append(products_)
                        reductions['error'].append(err_pot)
                        reductions['coeffs'].append(coeffs)
        
    if (verbosity):
            print(' End stability window---- \n')
    
    print('CLEANING RESULTS to avoid duplicates')
    
    cleaned_oxidations={}
    cleaned_oxidations['pot']=[]
    cleaned_oxidations['products']=[]
    cleaned_oxidations['error']=[]
    cleaned_oxidations['coeffs']=[]
    
    cleaned_reductions={}
    cleaned_reductions['pot']=[]
    cleaned_reductions['products']=[]
    cleaned_reductions['error']=[]
    cleaned_reductions['coeffs']=[]
    
    bucket=[]
    for iox in range(len(oxidations['pot'])): 
        found=False
        for ibucket in bucket :
            if ( ( np.abs(oxidations['pot'][iox]-ibucket[0])<1.e-5) and 
                 ( np.abs(oxidations['error'][iox]-ibucket[1])<1.e-5) ) :
                found=True
                break
        if found :
            pass
        else:
            bucket.append((oxidations['pot'][iox],oxidations['error'][iox]))
            cleaned_oxidations['pot'].append(oxidations['pot'][iox])
            cleaned_oxidations['products'].append(oxidations['products'][iox])
            cleaned_oxidations['error'].append(oxidations['error'][iox])
            cleaned_oxidations['coeffs'].append(oxidations['coeffs'][iox])
            
    bucket=[]
    for ired in range(len(reductions['pot'])): 
        found=False
        for ibucket in bucket :
            if ( ( np.abs(reductions['pot'][ired]-ibucket[0])<1.e-5) and 
                 ( np.abs(reductions['error'][ired]-ibucket[1])<1.e-5) ) :
                found=True
                break
        if found :
            pass
        else:
            bucket.append((reductions['pot'][ired],reductions['error'][ired]))
            cleaned_reductions['pot'].append(reductions['pot'][ired])
            cleaned_reductions['products'].append(reductions['products'][ired])
            cleaned_reductions['error'].append(reductions['error'][ired])
            cleaned_reductions['coeffs'].append(reductions['coeffs'][ired])        
        
    return cleaned_oxidations, cleaned_reductions
