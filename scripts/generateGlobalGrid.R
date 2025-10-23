library(dggridR)


res = 10

grid = dgconstruct(
  projection = 'ISEA',
  aperture = 3,
  topology = 'HEXAGON',
  res = res
)

# res = dg_closest_res_to_spacing(grid, spacing = 100, round = 'down', metric = T)
grid = dgsetres(grid, res)
dgearthgrid(grid, savegrid = sprintf('../data/global_grid_res%s.shp', res))

global.cell = data.frame(cell = getSpPPolygonsIDSlots(global), row.names = getSpPPolygonsIDSlots(global))
global = SpatialPolygonsDataFrame(global, global.cell)
writeOGR(global, "../data/global.shp", "", "ESRI Shapefile")