import ase.db
import pickle
from itertools import product

print("Collecting type III heterojunctions...")

def candidate(left, right):
    # Apply attr to each combination of band energies with given suffix
    def calculation(attr, suffix): return (
        attr(row,band)
            for row in (left,right)
            for band in ('vbm'+suffix,'cbm'+suffix)
    )

    # Check that a is in ascending order
    def ascending(a): return all(
        a[i]<a[i+1] for i in range(len(a)-1)
    )

    # Find the x - y*n closest to 0
    def fit(l, r):
        y = min(abs(l),abs(r))
        if y==0:
            return 0
        x = max(abs(l),abs(r))
        return ((x + y/2)%y) - y/2

    def fit_lat(l, r):
        (l_lat, r_lat) = (
            row.data['results-asr.structureinfo.json']['kwargs']['data']['spglib_dataset']['std_lattice']
                for row in (l, r)
        )
        return sum((
            abs(fit(l_lat[y][x], r_lat[y][x]))
                for x, y in product((0,1,2), repeat=2)
        ))


    # Only compare bands calculated with the same functional, prefer GW > HSE > PBE
    for (functional,suffix) in (('GW','_gw'),('HSE','_hse'),('PBE','')):
        if all(calculation(hasattr,suffix)):
            lvbm,lcbm,rvbm,rcbm = calculation(getattr,suffix)
            # Center bands on evac=0
            lvbm,lcbm,rvbm,rcbm = (
                lvbm - left.evac,
                lcbm - left.evac,
                rvbm - right.evac,
                rcbm - right.evac
            )
            if all((
                # Only accept thermodynamically stable materials
                left.thermodynamic_stability_level == 3 and right.thermodynamic_stability_level == 3, 
                left.is_magnetic == 0 and right.is_magnetic == 0, # Do not accept magnetic materials
                functional == 'HSE', # Only accept HSE functionals
                getattr(left,'class',None) == getattr(right,'class',None), # Materials must have the same phase
                left.spgnum == right.spgnum, # Materials must have the same space group
                left.crystal_type == right.crystal_type, # Materials must have the same crystal structure
                not ('C' in left.crystal_type), # Only accept A or AB type crystal structures 
                ascending((lvbm,lcbm,rvbm,rcbm)) # Gap is broken and bands are rising 
            )):
                return (
                    left.formula,                   #0
                    left.uid,                       #1
                    lvbm,                           #2
                    lcbm,                           #3
                    right.formula,                  #4
                    right.uid,                      #5
                    rvbm,                           #6
                    rcbm,                           #7
                    functional,                     #8
                    getattr(left,'class',None),     #9
                    left.spacegroup,                #10
                    left.crystal_type,              #11
                    rcbm - lvbm,                    #12
                    left.cell_area/right.cell_area, #13
                    fit_lat(left, right)             #14
                )

headers = [     # The title of each candidate column
    'Left',     #0
    'uid',      #1
    'vbm',      #2
    'cbm',      #3
    'Right',    #4
    'uid',      #5
    'vbm',      #6
    'cbm',      #7
    'Calc',     #8
    'Class',    #9
    'Group',    #10
    'Type',     #11
    'ΔE',       #12
    'R',        #13
    'ΔL'        #14
]

try:
    with open('typeIII.pickle','rb') as pk:
        typeIII = pickle.load(pk)
        print("Loaded from pickle.")
except (OSError, IOError):
    # Connect to the database
    db = ase.db.connect('c2db.db')
    # Find rows that have all of the parameters used to choose candidates
    rows = list(db.select('evac,thermodynamic_stability_level,is_magnetic,spgnum,crystal_type,has_asr_structureinfo'))
    # Find pairs of rows that could be type III heterojunctions
    typeIII = sorted(
        list( filter(
            lambda row: row != None, # Remove results that were not candidates
            [candidate(l,r) for l in rows for r in rows] # Take every pair in the selection
        ) ),
        key = lambda t: t[14] # Sort by how well the lattices fit
    )
    print("Loaded from ase db.")
    with open('typeIII.pickle','wb') as pk: # Cache results to make plot changes load faster
        pickle.dump(typeIII, pk)
        print("Dumped to pickle.")
