def sum_S_yield(community_members_array):
    sum_S_yield = ""
    S_yield = "- ( mu_#CELL# * N_#CELL# / $g#CELL#$ )"

    for strain in community_members_array:
        print(strain.name)
        if strain.consumes_S:
            sum_S_yield = sum_S_yield + S_yield.replace('#CELL#', strain.name)

    return sum_S_yield
