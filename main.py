from OCC.Extend.DataExchange import read_step_file
from part import Part
import argparse
import calculator

# CONSTANTS
# conversions for dimensional analysis
G_PER_KG = 1000
KG_PER_G = 0.001
MM3_PER_CM3 = 1000
CM3_PER_MM3 = 0.001
S_PER_HR = 3600
HR_PER_S = 1 / 3600

ALUMINUM_DENSITY_G_PER_CM3 = 2.7  # g/cm^3
ALUMINUM_PRICE_USD_PER_KG = 10.00  # $/kg
ALUMINUM_DENSITY_G_PER_MM3 = ALUMINUM_DENSITY_G_PER_CM3 * CM3_PER_MM3  # g/mm^3

COARSE_MILLING_MM3_PER_S = 500
FINE_MILLING_MM3_PER_S = 10

CNC_RATE_USD_PER_HR = 100.00
CNC_RATE_USD_PER_S = CNC_RATE_USD_PER_HR * HR_PER_S

CUT_FEE_USD = 10.00
SETUP_TIME_HR = 1.5
STOCK_FACTOR = 1.6  # rough estimate for typical 3d milling with fixturing and grip
MIN_BUY_MASS_G = 1000

# COMMAND LINE ARGS
# Absolute or relative path to the STEP file
file_path = ""
msg = "Pass a .step file as a command line argument to this program, analyzes and estimates the cost to manufacture part"
parser = argparse.ArgumentParser(description=msg)
parser.add_argument("filename", help="filepath to .step file to be analyzed")
args = parser.parse_args()

if args.filename:
    file_path = args.filename

# Make sure the file actually exists at this location
try:
    # Read the STEP file into a TopoDS_Shape
    my_shape = read_step_file(file_path)
except Exception as error:
    print(error)
    quit()

part = Part(my_shape)

# bounding box
P, V1, V2, V3 = part.getBoundingBox()
print(f"Bounding box dimensions: {V1.Magnitude():,.3f} mm × {V2.Magnitude():,.3f} mm × {V3.Magnitude():,.3f} mm")
bounding_box_volume_mm3 = V1.Magnitude() * V2.Magnitude() * V3.Magnitude()
# print(f"Bounding box volume: {bounding_box_volume_mm3:,.2f} mm^3")

# material cost = volume * density * price = mass * price
bounding_box_mass_g = bounding_box_volume_mm3 * ALUMINUM_DENSITY_G_PER_MM3
block_mass_kg = max(bounding_box_mass_g * STOCK_FACTOR, MIN_BUY_MASS_G) * KG_PER_G
# print(f"Block mass: {block_mass_kg:,.2f} kg")
block_volume_mm3 = block_mass_kg * G_PER_KG * (1 / ALUMINUM_DENSITY_G_PER_MM3)
block_material_cost_usd = block_mass_kg * ALUMINUM_PRICE_USD_PER_KG
# print(f"Block volume: {block_volume_mm3:,.2f} mm^3")
print(f"Material cost (block): ${block_material_cost_usd:,.2f}\n")

# volume
part_volume_mm3 = part.getVolumeMM3()  # mm^3
print(f"Part volume: {part_volume_mm3:,.2f} mm^3")

# surface area
part_surface_area = part.getSurfaceAreaMM2()
print(f"Surface Area: {part_surface_area:,.2f} mm^2")

# convex hull volume
convex_hull_volume_mm3 = part.getConvexHullVolumeMM3()
# print(f"Convex hull volume: {convex_hull_volume_mm3:,.2f} mm^3")

# milling volumes
coarse_volume_mm3 = block_volume_mm3 - convex_hull_volume_mm3
fine_volume_mm3 = convex_hull_volume_mm3 - part_volume_mm3

# print(f"Coarse volume : {coarse_volume_mm3:,.2f} mm^3")
# print(f"Fine volume : {fine_volume_mm3:,.2f} mm^3\n")

# milling times
coarse_milling_time_s = coarse_volume_mm3 / COARSE_MILLING_MM3_PER_S
fine_milling_time_s = fine_volume_mm3 / FINE_MILLING_MM3_PER_S
total_milling_time_s = coarse_milling_time_s + fine_milling_time_s
total_milling_time_hr = total_milling_time_s * HR_PER_S
total_milling_cost_usd = total_milling_time_s * CNC_RATE_USD_PER_S

# print(f"Coarse milling time : {coarse_milling_time_s:,.2f} s")
# print(f"Fine milling time : {fine_milling_time_s:,.2f} s")
print(f"Total milling time : {total_milling_time_hr:,.2f} hr")
# print(f"Total milling cost : ${total_milling_cost_usd:,.2f}\n")

setup_cost_usd = SETUP_TIME_HR * CNC_RATE_USD_PER_HR
total_cost = CUT_FEE_USD + setup_cost_usd + block_material_cost_usd + total_milling_cost_usd

print(f"TOTAL COST (1 ct): ${total_cost:,.2f}\n")

print("Qty  UnitPrice  Subtotal")
print(f"{1}    ${total_cost:,.2f}     ${total_cost:,.2f}")
for i in range(2, 6):
    unit_price = calculator.getUnitPriceUSD(i, total_cost)
    print(f"{i}    ${unit_price:,.2f}     ${unit_price * i:,.2f}")
for i in range(10, 26, 5):
    unit_price = calculator.getUnitPriceUSD(i, total_cost)
    print(f"{i}   ${unit_price:,.2f}      ${unit_price * i:,.2f}")

print(f"Complexity (0-1): {part.getComplexity():.3f}")
print(f"Lead Time: {calculator.getLeadTimeDays(part)} days")

