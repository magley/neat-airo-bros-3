-- https://datacrystal.romhacking.net/wiki/Super_Mario_Bros._3:RAM_map

local SAVESTATE_FNAME = "./SMB3_NeatAiro.State";
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

local INPUT_RADIUS = 7;
local INPUT_COUNT = ((INPUT_RADIUS*2 + 1) * (INPUT_RADIUS*2 + 1)) + 1;

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

-- Return array of integers each representing a value in the input grid.
-- 0 is empty space, 1 is solid block, -1 is enemy. Enemies override solids.
local function get_inputs()
    -- We fill evereything with 0 -> empty cell
    -- For each cell, check if it's a tile, and if so, give it a value of 1

    local inputs = {};

    for cy = -INPUT_RADIUS, INPUT_RADIUS do
        for cx = -INPUT_RADIUS, INPUT_RADIUS do
            local dy = cy * 16;
            local dx = cx * 16;

            inputs[#inputs + 1] = 0; -- Default cell value is 0 => nothing.

            local tile = get_tile(dx, dy);
            if tile == 1 then
                inputs[#inputs] = 1; -- 1 => solid block/floor/wall.
            end
        end
    end
    return inputs;
end

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
            gui.drawBox(
                cell.x - RENDER_NODE_SIZE/2,
                cell.y - RENDER_NODE_SIZE/2,
                cell.x + RENDER_NODE_SIZE/2,
                cell.y + RENDER_NODE_SIZE/2,
                0xFF000000,
                0xFFFFFFFF
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

    emu.frameadvance();
end