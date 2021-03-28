from tabulate import tabulate
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import ase.db

print("Collecting type III heterojunctions from c2db.db...")

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

    # Only compare bands calculated with the same functional, prefer GW > HSE > PBE
    for (functional,suffix) in (('GW','_gw'),('HSE','_hse'),('PBE','')):
        if all(calculation(hasattr,suffix)):
            lvbm,lcbm,rvbm,rcbm = calculation(getattr,suffix)
            lvbm,lcbm,rvbm,rcbm = ( # Center bands on evac=0
                lvbm - left.evac,
                lcbm - left.evac,
                rvbm - right.evac,
                rcbm - right.evac
            )
            if all((
                functional != 'PBE', # Do not accept PBE functionals
                getattr(left,'class',None) == getattr(right,'class',None), # Materials must have the same phase
                left.spgnum == right.spgnum, # Materials must have the same space group
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
                    rcbm - lvbm,                    #11
                    left.cell_area/right.cell_area  #12
                )

try:
    with open('./typeIII.pickle','rb') as pk:
        typeIII = pickle.load(pk)
        print("Loaded from pickle.")
except (OSError, IOError):
    # Connect to database
    db = ase.db.connect('c2db.db')
    # Find rows that have a vacuum energy, space group number, high thermodynamic stability, and are not magnetic
    rows = list(db.select('evac,spgnum,thermodynamic_stability_level=3,is_magnetic=0'))
    # Find pairs of rows that could be type III heterojunctions
    typeIII = sorted(
        filter(
            None.__ne__, # Remove results that were not candidates
            [candidate(l,r) for l in rows for r in rows] # Take every pair in the selection
        ),
        key = lambda t: t[12] if 1 <= t[12] else 1/t[12] # Sort by how close R is to 1
    )
    print("Loaded from ase db.")
    with open('./typeIII.pickle','wb') as pk: # Cache results to make plot changes load faster
        pickle.dump(typeIII, pk)
        print("Dumped to pickle.")

spgs = pd.unique([t[10] for t in typeIII]) # Find which space groups we found pairs from
print("Collected " + str(len(typeIII)) + " type III heterojunctions from " + str(len(spgs)) + " space groups:")
print(spgs)

best = [next(filter(lambda t: t[10] == spg, typeIII)) for spg in spgs] # Only plot the best pair from each spacegroup
headers = [
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
    'ΔE',       #11
    'R'         #12
]
table_headers_indices = [0,4,8,10,11,12] # Which headers to display in the table
print(tabulate(best, headers)) # Print a table of the pairs with their UIDs

def bands_with_table(pairs): # Plot each pair of candidates, with bars for the band gaps, and a table for other stats
    ax = plt.subplot(label=str(id(pairs)))
    ax.margins(x=0)

    flatx = range(len(pairs)*2)

    leftvbm = np.array([t[2] for t in pairs])
    leftcbm = np.array([t[3] for t in pairs])
    emin = np.amin(leftvbm) - .5

    rightvbm = np.array([t[6] for t in pairs])
    rightcbm = np.array([t[7] for t in pairs])
    emax = np.amax(rightcbm) + .5

    leftx = flatx[::2]
    ax.bar(leftx, leftvbm-emin, bottom=emin, color='C0')
    ax.bar(leftx, leftcbm-emax, bottom=emax, color='C1')
    #ax.bar(leftx, leftvbm-leftcbm, bottom=leftcbm)

    rightx = flatx[1::2]
    ax.bar(rightx, rightvbm-emin, bottom=emin)
    ax.bar(rightx, rightcbm-emax, bottom=emax)
    #ax.bar(rightx, rightvbm-rightcbm, bottom=rightcbm)

    ax.set_ylim(emin,emax)

    plt.xticks([],[])

    table = plt.table(
        cellText = [
            [ t[i] if not isinstance(t[i], float) else round(t[i],5) for t in pairs ]
                for i in table_headers_indices
        ],
        rowLabels=[headers[i] for i in table_headers_indices],
        loc='bottom',
        picker=True # Table is clickable
    )
    table.auto_set_font_size(False)
    table.set_fontsize('small')

    plt.axhline(0, linewidth=.75, color='k')

    for v in leftx:
        plt.axvline(v-.5, linewidth=1.25, color='k')

    plt.ylabel('Band gap, Blue=VBM, Orange=CBM (eV)')

fig = plt.figure(figsize=(15,7))
bands_with_table(best)
plt.title("Top candidates (" + str(len(spgs)) + " space groups)")
plt.tight_layout()

def pick_handler(event): # Table was clicked
    (tx0, ty0, tw, th) = event.artist.get_window_extent(fig.canvas.renderer).bounds
    (cx0, cy0, cw, ch) = event.artist[len(table_headers_indices)-1,-1].get_window_extent(fig.canvas.renderer).bounds
    x = len(best)*(event.mouseevent.x-cx0-cw)/(tw-cw)
    y = len(table_headers_indices)*(event.mouseevent.y-cy0)/th

    spg_index = int(x) # Which column was clicked
    spg_best = [t for t in typeIII if t[10] == spgs[spg_index]] # Rows in this group
    top_index = 0
    top = spg_best[top_index:top_index+20] # Top 20 rows in this group

    if(y < 2): # Scatter plot of Delta E and R
        plt.figure()
        plt.scatter(
            x=[s[11] for s in top],
            y=[s[12] for s in top]
        )
        for s in top: # Label scatter plot points with formulas
            plt.annotate(s[0]+","+s[4],(s[11],s[12]))
        plt.xlabel(headers[11] + " (eV)")
        plt.ylabel(headers[12] + " (1)")
    else: # VBM,CBM bars for this group
        plt.figure(figsize=(15,7))
        bands_with_table(top)

    plt.title(spgs[spg_index] + " (top " + str(len(top)) + " of " + str(len(spg_best)) + " candidates)")
    plt.tight_layout()
    plt.get_current_fig_manager().toolbar.children['!checkbutton'].pack_forget()
    plt.get_current_fig_manager().toolbar.children['!checkbutton2'].pack_forget()
    plt.show(block=False)

plt.get_current_fig_manager().toolbar.children['!checkbutton'].pack_forget()
plt.get_current_fig_manager().toolbar.children['!checkbutton2'].pack_forget()
fig.canvas.mpl_connect('pick_event',pick_handler)

fig.savefig('fig1.png')
plt.show()
