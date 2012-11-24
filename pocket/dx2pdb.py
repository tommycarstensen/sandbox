counts = lines[0].strip().split()[-3:]
origin = lines[1].split()[-3:]
for i in range(3):
    counts[i] = int(counts[i])
    origin[i] = float(origin[i])
delta_x = float(lines[2].split()[1])
delta_y = float(lines[3].split()[2])
delta_z = float(lines[4].split()[3])
delta = [delta_x,delta_y,delta_z,]

for line in lines:

    x_grid = ( i%(counts[2]*counts[1]*counts[0]) - i%(counts[2]*counts[1]) ) / (counts[2]*counts[1])
    y_grid = ( i%(counts[2]*counts[1]          ) - i% counts[2]            ) /  counts[2]
    z_grid =   i% counts[2]
    x = origin[0] + delta[0]*x_grid
    y = origin[1] + delta[1]*y_grid
    z = origin[2] + delta[2]*z_grid
