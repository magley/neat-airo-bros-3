-------------------------------------------------------------------------------
--
-- NEAT AIRO BROS 3
-- Neuro Evolution of Augmenting Topologies for Super Mario Bros. 3
--
-------------------------------------------------------------------------------


-- https://datacrystal.romhacking.net/wiki/Super_Mario_Bros._3:RAM_map

local SAVESTATE_FNAME = "./SMB3_NeatAiro_lvl3.State";
local TILE_BLOCK_VALUES = {
    0xe2, 0xe3, 0xe4, -- Blue jumpthru block
    0xa0, 0xa0, 0xa2, -- Green jumpthru block
    0x25, 0x26, 0x27, -- White jumpthru block
    0x61, 0x62, 0x63, -- ? Blocks
    0x50, 0x51, 0x52, -- Red jumpthru block
    0x53, 0x54, 0x55, -- Wood top
    0x56, 0x57, 0x58, -- Wood middle
    0xad, 0xae,
    0xb1, 0xb2, -- Pipe
    0xba, 0xbb, -- Pipe
    0xa0, 0xa1, 0xa2,
    0x57,
    0x67,
    0x79,
    0x2c,
    0x5f,
    -- Special blocks (note blocks, item bricks, movable wood)
    103, 109, 107, 106, 105, 46, 116, 115,
};

local ENEMY_VALUES = {
    0x6C, -- green koopa troopa
    0x6D, -- red koopa troopa
    0x6E, -- green koopa paratroopa jumpy
    0x72, -- goomba
    0x73, -- para goomba (red)
    0xA0, -- green piranha plant up
    0x82, -- boomerang bro
};

local BUTTONS = {
    "A",
    "B",
    "Left",
    "Right",
}

local INPUT_RADIUS = 7;
local INPUT_COUNT = ((INPUT_RADIUS*2 + 1) * (INPUT_RADIUS*2 + 1)) + 1;
local OUTPUT_COUNT = #BUTTONS;

local MAX_NODES = 1000; -- INPUT_COUNT + OUTPUT_COUNT + hidden count

local tick = 0;
local mario_x = 0;
local mario_y = 0;

-------------------------------------------------------------------------------
-- Mario 3 RAM-related functions (network input layer).
-------------------------------------------------------------------------------

-- Update `mario_x` and `mario_y` based on the player's location in-game.
local function update_player_position()
    -- SMB3 uses pages, a page is a chunk of 256 pixels.
    -- position = (current_page, position_in_current_page).

    local x_in_page = memory.readbyte(0x90); -- x % 256
    local y_in_page = memory.readbyte(0xA2); -- y % 256

    local x_page = memory.readbyte(0x75); -- floor(x / 256)
    local y_page = memory.readbyte(0x87); -- floor(y / 256)

    local x = x_in_page + x_page * 256;
    local y = y_in_page + y_page * 256;

    mario_x = x
    mario_y = y
end

-- Return byte value in tilemap at cell located at Mario's position + (dx, dy).
local function get_tile_raw(dx, dy)
    -- The level's tilemap begins at 0x6000.
    -- Each tile is a 1 byte value.
    -- Tiles are ordered in 16 block wide rows.
    -- Each page has 27 rows (applies to 1-1, 1-3, not sure about the others).

    local BASE_ADDR = 0x6000;
    local COLS_IN_PAGE = 16;
    local ROWS_IN_PAGE = 27;
    local CELLS_IN_PAGE = COLS_IN_PAGE * ROWS_IN_PAGE;

    local x = mario_x + dx + 8;  -- +8 because mario's position is top-left.
    local y = mario_y + dy;

    local cell_x = math.floor(x / 16);
    local cell_y = math.floor(y / 16);

    local row = cell_y % ROWS_IN_PAGE;
    local col = cell_x % COLS_IN_PAGE;
    local page = math.floor(cell_x / COLS_IN_PAGE);

    local addr = BASE_ADDR + (CELLS_IN_PAGE * page) + (row * COLS_IN_PAGE) + (col);
    local tile = memory.readbyte(addr);

    return tile
end

-- Return `true` if `raw_byte_value` is the value of a solid block.
-- Solid blocks are defined in `TILE_BLOCK_VALUES`.
local function is_solid_block(raw_byte_value)
    for i=1, #TILE_BLOCK_VALUES do
        if TILE_BLOCK_VALUES[i] == raw_byte_value then
            return true
        end
    end
    return false
end

-- `1` if tile at cell located at Mario's position + (dx, dy) is solid, else `0`.
local function get_tile(dx, dy)
    local tile_value_raw = get_tile_raw(dx, dy)
    if is_solid_block(tile_value_raw) then
        return 1
    else
        return 0
    end
end

-- Return `true` if `raw_byte_value` is the value of a an enemy.
-- Enemies are defined in `ENEMY_VALUES`.
local function is_enemy(raw_byte_value)
    for i=1, #ENEMY_VALUES do
        if ENEMY_VALUES[i] == raw_byte_value then
            return true
        end
    end
    return false
end

-- Return array of objects `{x, y, id}`.
-- Each object is a cell denoting a single enemy.
local function get_enemies()
    -- Enemy data begins at 0x7B40.
    -- The first byte is always 1.
    -- 
    -- Each enemy has 3 bytes associated to it:
    -- EnemyID, StartX, StartY.
    -- EnemyID determines the type of enemy.
    -- StartX and StartY are in tile cells.
    -- 
    -- After the last enemy on screen, a 0xFF byte follows
    -- and then 0x00 bytes until the end of enemy buffer.
    -- 
    -- The positions are split into 2 values, much like
    -- Mario's position. The game supports 5 enemies at
    -- once (???) and their positions go in reverse order
    -- in the RAM map. The game recycles enemy locations
    -- by [i mod 5].
    -- 
    -- Enemy state is stored in 0x661 - 0x665.
    -- if 0 => dead
    -- if 6 => knocked by shell (flies down the screen).
    -- if 7 => stopmed
    -- if 8 => killed by shell (turns to dust).

    local BASE_ENEMY_ID = 0x7B40 + 0x01;
    local BASE_X_HIGH = 0x76;
    local BASE_X_LOW = 0x91;
    local BASE_Y_HIGH = 0x88;
    local BASE_Y_LOW = 0xA3;
    local BASE_ENEMY_STATE = 0x661;
    local enemies = {};

    for i = 0, 4 do
        local addr = BASE_ENEMY_ID + i * 3;

        local enemy_type = memory.readbyte(addr);
        local x_in_page = memory.readbyte(BASE_X_LOW + (4-i) % 5);
        local y_in_page = memory.readbyte(BASE_Y_LOW + (4-i) % 5);
        local x_page = memory.readbyte(BASE_X_HIGH + (4-i) % 5);
        local y_page = memory.readbyte(BASE_Y_HIGH + (4-i) % 5);
        local state = memory.readbyte(BASE_ENEMY_STATE + (4-i) % 5);

        local x = x_in_page + x_page * 256;
        local y = y_in_page + y_page * 256;

        if is_enemy(enemy_type) and (state < 6 and state > 0) then
            if enemy_type == ENEMY_VALUES[7] then
                y = y + 8;
            end

            enemies[#enemies + 1] = {
                ["id"] = enemy_type,
                ["x"] = x,
                ["y"] = y,
            };
        end
    end

    return enemies;
end


-- Return array of integers each representing a value in the input grid.
-- 0 is empty space, 1 is solid block, -1 is enemy. Enemies override solids.
-- **Does not** include the bias node.
local function get_inputs()
    -- We fill evereything with 0 -> empty cell
    -- For each cell, check if it's a tile, and if so, give it a value of 1.
    -- For each enemy, check if its on the i-th cell and if so, give it -1.

    local inputs = {};

    local enemies = get_enemies();

    for cy = -INPUT_RADIUS, INPUT_RADIUS do
        for cx = -INPUT_RADIUS, INPUT_RADIUS do
            local dy = cy * 16;
            local dx = cx * 16;

            inputs[#inputs + 1] = 0; -- Default cell value is 0 => nothing.

            for i = 1, #enemies do
                local enem_x = enemies[i].x - 8;
                local enem_y = enemies[i].y;

                local dist_x = math.abs(enem_x - (mario_x + dx));
                local dist_y = math.abs(enem_y - (mario_y + dy));

                if dist_x <= 8 and dist_y <= 8 then
                    inputs[#inputs] = -1;
                end
            end

            local tile = get_tile(dx, dy);
            if tile == 1 then
                inputs[#inputs] = 1; -- 1 => solid block/floor/wall.
            end
        end
    end
    return inputs;
end

-------------------------------------------------------------------------------
-- Math
-------------------------------------------------------------------------------

-- For some reason, `math.tanh` was deprecated.
-- http://lua-users.org/wiki/HyperbolicFunctions
local function tanh(x)
    if x == 0 then return 0.0 end
    local neg = false
    if x < 0 then x = -x; neg = true end
    if x < 0.54930614433405 then
      local y = x * x
      x = x + x * y *
          ((-0.96437492777225469787e0  * y +
            -0.99225929672236083313e2) * y +
            -0.16134119023996228053e4) /
          (((0.10000000000000000000e1  * y +
             0.11274474380534949335e3) * y +
             0.22337720718962312926e4) * y +
             0.48402357071988688686e4)
    else
      x = math.exp(x)
      x = 1.0 - 2.0 / (x * x + 1.0)
    end
    if neg then x = -x end
    return x
  end

local function sigmoid(x)
    return tanh(x)
end

-------------------------------------------------------------------------------
-- Data structures
-------------------------------------------------------------------------------

-- A MutationRate is a mutable record of probabilities for mutation.
-- 1 genome <--> 1 mutation rate.
local function MutationRate()
    local mutation_rate = {};
    mutation_rate.connections = 0.25;
    mutation_rate.link = 2.0;
    mutation_rate.bias = 0.4;
    mutation_rate.node = 0.5;
    mutation_rate.enable = 0.2;
    mutation_rate.disable = 0.4;
    mutation_rate.step = 0.1;

    return mutation_rate;
end

-- Single node in a NN.
local function Neuron()
    local neuron = {};
    neuron.incoming = {}; -- Incoming genes.
    neuron.value = 0;

    return neuron;
end

-- Connection between 2 neurons.
local function Gene()
    local gene = {};
    gene.into = 0;
    gene.out = 0;
    gene.weight = 0;
    gene.enabled = true;
    gene.innovation = 0;

    return gene;
end

-- Create deep copy of `gene`.
local function GeneCopy(gene)
    local new = Gene()
	new.into = gene.into
	new.out = gene.out
	new.weight = gene.weight
	new.enabled = gene.enabled
	new.innovation = gene.innovation

    return new;
end

-- A Genome is an individual competing in a population.
local function Genome()
    local genome = {};
    genome.genes = {};
    genome.fitness = 0;
    genome.fitness_adj = 0;
    genome.network = {};
    genome.max_neuron = 0;
    genome.global_rank = 0;
    genome.mutation_rate = MutationRate();

    return genome;
end

-- Create deep copy of `genome`.
local function GenomeCopy(genome)
    local new = Genome();

    for i = 1, #genome.genes do
        table.insert(new.genes, GeneCopy(genome.genes[i]));
    end

    for mut, rate in pairs(genome.mutation_rate) do
        new.mutation_rate[mut] = rate;
    end

    new.genome = genome.max_neuron;

    return new;
end

-- A Species is a collection of genomes grouped by similar features.
local function Species()
    local species = {};
    species.genomes = {};
    species.staleness = 0;
    species.top_fitness = 0;
    species.avg_fitness = 0;

    return species;
end

-- A Pool is a collection of species spanning multiple generations.
local function Pool()
    local pool = {};
    pool.species = {};
    pool.generation = 0;
    pool.innovation = OUTPUT_COUNT;
    pool.curr_species = 1;
    pool.curr_genome = 1;
    pool.top_fitness = 0;

    return pool;
end
local pool = {};

-- Modifies the global pool object.
local function new_innovation()
    pool.innovation = pool.innovation + 1;
    return pool.innovation;
end

-------------------------------------------------------------------------------
-- Neural network high level functions.
-------------------------------------------------------------------------------

-- Create a network from the given `genome: Genome`.
-- Returns the created network.
-- You probably want `generate_network_for()` instead.
local function Network(genome)
    local network = {};
    network.neurons = {};

    -- Initialize all input and output neurons.

    for i = 1, INPUT_COUNT do
        network.neurons[i] = Neuron();
    end
    for j = 1, OUTPUT_COUNT do
        network.neurons[MAX_NODES + j] = Neuron();
    end

    -- Sort by `out` i.e. semi-horizontal ordering.

    table.sort(genome.genes,
        function(g1, g2)
            return (g1.out < g2.out)
        end
    );

    -- For each enabled gene in the genome, add its neurons and connect them.

    for i = 1, #genome.genes do
        local g = genome.genes[i];

        if g.enabled then
            if network.neurons[g.out] == nil then
                network.neurons[g.out] = Neuron();
            end
            table.insert(network.neurons[g.out].incoming, g);
            if network.neurons[g.into] == nil then
                network.neurons[g.into] = Neuron();
            end
        end
    end

    return network;
end

-- Given the genome, create its network.
-- Modifies genome.
local function generate_network_for(genome)
    local network = Network(genome);
    genome.network = network;
end

-- Feed `input_data` to `network` and spit out `output`.
-- Returns a map between button names and true/false, for each button.
--
-- Modifies `network`: updates its input nodes to `input_data`.
--
-- Modifies `input_data`: adds a bias node.
local function network_evaluate(network, input_data)
    -- Add bias node.
    table.insert(input_data, 1);

    -- Just in case.
    if #input_data ~= INPUT_COUNT then
        console.writeline("Input count expected " .. INPUT_COUNT .. " got " .. #input_data);
        return {};
    end

    -- Evaluate input layer.
    for i = 1, INPUT_COUNT do
        network.neurons[i].value = input_data[i];
    end

    -- Evaluate output and hidden layer.
    for _, neuron in pairs(network.neurons) do -- `neurons` is sparse
        local sum = 0;

        for i = 1, #neuron.incoming do
            local gene = neuron.incoming[i];
            local neuron_in = network.neurons[gene.into];

            sum = sum + gene.weight * neuron_in.value;
        end

        if #neuron.incoming > 0 then -- Do this only for non-inputs.
            neuron.value = sigmoid(sum);
        end
    end

    -- Evaluate game action based on network state.
    local output_arr = {}; -- true/false for each button.
    for i = 1, OUTPUT_COUNT do
        local btn_name = "P1 " .. BUTTONS[i];

        if network.neurons[MAX_NODES + i].value > 0 then
            output_arr[btn_name] = true;
        else
            output_arr[btn_name] = false;
        end
    end

    return output_arr;
end

-------------------------------------------------------------------------------
-- Mutation
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
-- Render
-------------------------------------------------------------------------------

local RENDER_NODE_SIZE = 5;
local RENDER_LEFT = 10;
local RENDER_TOP = 10;

-- These are centered values.

local RENDER_X = RENDER_LEFT + (INPUT_RADIUS) * RENDER_NODE_SIZE - RENDER_NODE_SIZE/2;
local RENDER_Y = RENDER_TOP + (INPUT_RADIUS) * RENDER_NODE_SIZE - RENDER_NODE_SIZE/2;

-- Returns dictionary of cells
-- `{x, y, value}`
local function draw_get_cells()
    local cells = {};

    local inputs = get_inputs();

    -- Input nodes

    local k = 1;
    for dy = -INPUT_RADIUS, INPUT_RADIUS do
        for dx = -INPUT_RADIUS, INPUT_RADIUS do
            local cell = {};
            cell.x = RENDER_X + RENDER_NODE_SIZE * dx;
            cell.y = RENDER_Y + RENDER_NODE_SIZE * dy;
            cell.value = inputs[k];
            cells[k] = cell;

            k = k + 1;
        end
    end

    -- Bias node

    local bias_cell = {};
    bias_cell.x = RENDER_X + (INPUT_RADIUS) * RENDER_NODE_SIZE;
    bias_cell.y = RENDER_Y + (INPUT_RADIUS+1) * RENDER_NODE_SIZE;
    bias_cell.value = inputs[INPUT_COUNT];
    cells[INPUT_COUNT] = bias_cell;

    return cells;
end

local function draw()
    -- Input box

    gui.drawBox(
        RENDER_X - (INPUT_RADIUS) * RENDER_NODE_SIZE - RENDER_NODE_SIZE/2,
        RENDER_Y - (INPUT_RADIUS) * RENDER_NODE_SIZE - RENDER_NODE_SIZE/2,
        RENDER_X + (INPUT_RADIUS) * RENDER_NODE_SIZE + RENDER_NODE_SIZE/2,
        RENDER_Y + (INPUT_RADIUS) * RENDER_NODE_SIZE + RENDER_NODE_SIZE/2,
        0xFF00002F,
        0x708080CD
    );

    -- Input nodes

    local cells = draw_get_cells();
    for n, cell in pairs(cells) do
        -- For input nodes, draw only non-zero ones.
        -- For hidden/output nodes, always draw them.
        if n > INPUT_COUNT or cell.value ~= 0 then
            -- Determine color based on value (block/enemy)
            local color = 0xFFFFFFFF;
            if cell.value == -1 then
                color = 0xFF000000;
            end


            gui.drawBox(
                cell.x - RENDER_NODE_SIZE/2,
                cell.y - RENDER_NODE_SIZE/2,
                cell.x + RENDER_NODE_SIZE/2,
                cell.y + RENDER_NODE_SIZE/2,
                0xFF000000,
                color
            );
        end
    end

    -- Player

    gui.drawBox(
        RENDER_X - RENDER_NODE_SIZE/4,
        RENDER_Y + RENDER_NODE_SIZE - RENDER_NODE_SIZE/2,
        RENDER_X + RENDER_NODE_SIZE/4,
        RENDER_Y + 2 * RENDER_NODE_SIZE - RENDER_NODE_SIZE/2,
        0xFF660000,
        0xFFDD0000
    );
end

-------------------------------------------------------------------------------
-- Program
-------------------------------------------------------------------------------

local function on_create()
    savestate.load(SAVESTATE_FNAME);
end

local function on_update()
    update_player_position();
    get_inputs();

    draw();
end

on_create();
while true do
    on_update();

    tick = tick + 1;
    emu.frameadvance();
end
