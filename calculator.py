# rough reverse engineer from xometry
def getUnitPriceUSD(count, init_price):
    multiplier = 0.000244734 * pow(init_price, 2) + 0.788575 * init_price + 41.8959
    exp = -0.00000234429 * pow(init_price, 2) + 0.00327053 * init_price - 1.28692
    unit_price = multiplier * pow(count, exp)
    return unit_price

# rough reverse engineer from xometry
def getLeadTimeDays(part, base_days=3, k_volume=1/200000, k_dimension=1/150):
    volume_mm3 = part.getVolumeMM3()
    max_dim_mm = part.getMaxDimMM()
    lead_time = base_days + (k_volume * volume_mm3) + (k_dimension * max_dim_mm)
    # Round up to nearest whole day
    return int(round(lead_time))