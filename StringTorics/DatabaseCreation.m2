----------------------------------
-- Code for creating data bases --
----------------------------------

createPolytopeDatabase = method(
    Options => {
        "Hodge" => null, -- TODO: not used yet
        "Count" => null  -- TODO: not used yet
        })

createPolytopeDatabase(String, List) := opts -> (dbfilename, topes) -> (
    -- open data base file
    F := openDatabaseOut dbfilename;
    -- F#"info" = "4990 reflexive polytopes of h11=5"
    -- F#"topes" = toString topes;
    -- loop through topes, create CYPolytopeData, populate it, write it to data base.
    elapsedTime for i from 0 to #topes - 1 do elapsedTime (
        << "computing for polytope " << i << endl;
        V := cyPolytopeData(topes#i, ID => i); -- note that the polytope data is really that of the dual to topes#i.
        -- now fill it with data we want
        basisIndices V; -- compute them
        isFavorable V; -- compute h11, h21, favorability.
        annotatedFaces V; -- compute annotated faces
        -- now write it
        F#(toString i) = dump V;
        );
    --if opts#"Hodge" =!= null then F#"hodge" = toString(opts.Hodge);
    --if opts#"Count" =!= null then F#"count" = #topes;
    close F;
    )

addToCYDatabase = method(Options => {NTFE => false})

addToCYDatabase(String, CYPolytopeData) := opts -> (dbfilename, Q) -> (
    elapsedTime Xs := findAllCYs Q;
    << "  " << #Xs << " triangulations total" << endl;
    if opts.NTFE then (
        elapsedTime H := partition(restrictTriangulation, Xs);
        << "  " << #(keys H) << " NTFE triangulations" << endl;
        Xs = (keys H)/(k -> H#k#0); -- only take one triangulation that matches
        );
    F := openDatabaseOut dbfilename;
    for X in Xs do (
        F#(toString label X) = dump X;
        );
    close F;    
    )
addToCYDatabase(String, Database, ZZ) := opts -> (dbfilename, topesDB, i) -> (
    <<  "-- doing polytope " << i << endl;
    Q := cyPolytopeData(topesDB#(toString i), ID => i);
    addToCYDatabase(dbfilename, Q, opts);
    )
