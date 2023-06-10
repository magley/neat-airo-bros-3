local SAVESTATE_FNAME = "./SMB3_NeatAiro.State";


function OnCreate()
    savestate.load(SAVESTATE_FNAME);
end


OnCreate();