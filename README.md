# neat-airo-bros-3

NEAT Algorithm that plays Super Mario Bros. 3 for the NES.
Inspired by SethBling.

### Research info:

**Goal**: Train an agent to complete a single SMB3 level using NEAT.

**Result**: Unsuccessful, the agent gets stuck near the end.

**Analysis**: One of the production tests resulted in a success, but I was
not able to replicate it. This all boils down to RNG and parameter tuning.

For more info as well as a summary of NEAT, check the [research poster](https://github.com/magley/neat-airo-bros-3/blob/master/poster.pdf).

### How to run:

- open the BizHawk emulator
- open your legally obtained Super Mario Bros. 3 Rom file
- load the script: `Tools > Lua Console > Open Script > NeatAiro3.lua`
- create a savestate called `SMB3_NeatAiro_lvl3.State` in the same folder `NeatAiro3.lua` is in (do this just once)
