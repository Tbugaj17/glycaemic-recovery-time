-- PDF tables: (1) a line under headings that span several columns,
-- (2) tables kept on one page: a plain tabular instead of a longtable that
--     can break across pages (Quarto places it in its own table float with
--     the caption). Very long tables (> MAX_ROWS rows) are left as
--     longtables so they can still continue on the next page.
local MAX_ROWS = 45

local function cmidrules(row)
  local col, out = 1, {}
  for _, cell in ipairs(row.cells) do
    local span = cell.col_span or 1
    if (span > 1 or col > 1) and pandoc.utils.stringify(cell.contents) ~= "" then
      table.insert(out, string.format("\\cmidrule(lr){%d-%d}", col, col + span - 1))
    end
    col = col + span
  end
  return table.concat(out)
end

local function count_rows(tbl)
  local n = 0
  for _, body in ipairs(tbl.bodies) do n = n + #body.body end
  return n
end

function Table(tbl)
  if not FORMAT:match("latex") then return nil end
  local latex = pandoc.write(pandoc.Pandoc({tbl}), "latex")

  -- (1) cmidrules under spanning header cells
  local rows = tbl.head.rows
  if #rows >= 2 then
    local rules = cmidrules(rows[1])
    if rules ~= "" then
      local out, pos = {}, 1
      while true do
        local a, b = latex:find("\\toprule\\noalign{}", pos, true)
        if not a then break end
        local c, d = latex:find("\\\\\n", b, true)
        if not c then break end
        table.insert(out, latex:sub(pos, d)); table.insert(out, rules .. "\n")
        pos = d + 1
      end
      table.insert(out, latex:sub(pos)); latex = table.concat(out)
    end
  end

  -- (2) keep on one page: longtable -> tabular
  if count_rows(tbl) <= MAX_ROWS then
    local spec_start = latex:find("\\begin{longtable}[]{", 1, true)
    if spec_start then
      local open = spec_start + #"\\begin{longtable}[]"
      local depth, i = 0, open
      repeat
        local ch = latex:sub(i, i)
        if ch == "{" then depth = depth + 1 elseif ch == "}" then depth = depth - 1 end
        i = i + 1
      until depth == 0
      local spec = latex:sub(open, i - 1)
      local rest = latex:sub(i)
      -- caption (if any) is the first \caption{...} before \toprule
      local caption = ""
      local cs = rest:find("\\caption{", 1, true)
      local tr = rest:find("\\toprule", 1, true)
      if cs and tr and cs < tr then
        local d, j = 0, cs + #"\\caption"
        repeat
          local ch = rest:sub(j, j)
          if ch == "{" then d = d + 1 elseif ch == "}" then d = d - 1 end
          j = j + 1
        until d == 0
        caption = rest:sub(cs, j - 1)
        rest = rest:sub(tr)
      end
      -- keep only the first copy of the header, drop longtable markers
      local eh = rest:find("\\endfirsthead", 1, true)
      if eh then
        local h2 = rest:find("\\endhead", eh, true)
        rest = rest:sub(1, eh - 1) .. rest:sub(h2 + #"\\endhead")
      else
        rest = rest:gsub("\\endhead\n?", "", 1)
      end
      rest = rest:gsub("\\bottomrule\\noalign{}\n\\endlastfoot\n?", "", 1)
      rest = rest:gsub("\\end{longtable}", "\\bottomrule\\noalign{}\n\\end{tabular}", 1)
      latex = "\\begin{tabular}" .. spec .. "\n" .. rest
    end
  end
  return pandoc.RawBlock("latex", latex)
end
