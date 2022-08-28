# MC-NPC
Monte-Carlo Nuclear Pore Complex Simulation.

This code is implemented in a few functions, where the relationships between the different functions are shown in the pseudocode below. `GenerateLabelSites` and `GenerateFluorophoreSites` are implemented concurrently in the actual MATLAB code to capitalize on MATLAB's vectorization capabilities. Run `generateNPCs` for an example of the code.

# Pseudocode
```
OBTAIN Nup species and membrane size
CALL GetNupParameters with Nup species RETURNING location, radial distance, probability of labeling, and average SMLs per site
CALL GenerateCentroids with membrane size RETURNING x,y-coordinates of NPC centroids
FOR each NPC centroid
    CALL GenerateLabelSites with Nup species and location RETURNING list of x,y-coordinates of antibody label sites relative to the NPC centroid
    FOR each antibody label site
        CALL GenerateFluorophoreSite with location RETURNING list of x,y-coordinates of fluorophore site relative to the antibody label site
        CALL GenerateLocalizations with Nup species RETURNING list of x,y-coordinates of fluorophore localizations relative to the fluorophore location
        FOR each fluorophore blinking event
            COMPUTE absolute fluorophore x,y-coordinates by adding up the x,y-coordinates from above
        END FOR
    END FOR
END FOR
PRINT absolute fluorophore x,y-coordinates of single-molecule localizations for the whole membrane
```

# License
Distributed under the GNU General Public License v3.0. See `LICENSE.txt` for more information.

# Contact
Written by Wei Hong Yeo, [Functional Optical Imaging Laboratory](http://foil.northwestern.edu/), Northwestern University, 2022.