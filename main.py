# from tqdm import tqdm
import multiprocessing as mp
# import numpy as np
import math
import re

CFG_TARGET_SIMILARITY = 0.80 # [0.5 .. 1.0]
CFG_DATABASE_PART = 1
CFG_POOLS_CNT = 7

# database = {}
# database['human'] = 'bebrabebrabra'
# database['gigachad'] = 'bebrasupermegagigaberebebratheendbra'
# small_map = {}
# small_map['humanity'] = 'bebraabebra'
# small_map['coolness'] = 'giga'
def readChainsMap(path):
	chainMap = {}
	with open (path, "r") as file:
		content = file.read()
		titles = re.findall('>[^\n]*', content)
		chains = re.split('>[^\n]*\n', content)
		chains = chains[1:]
		for i, chain in enumerate(chains):
			chains[i] = chain.replace('\n', '')
		assert len(titles) == len(chains)
		for i in range(len(titles)):
			chainMap[titles[i]] = chains[i]
	return chainMap

# small_map = readChainsMap("lab0.txt")
database = readChainsMap("uniprot_sprot.fasta")
database_pools = []
for i in range(CFG_POOLS_CNT):
	database_pools.append({})
# chains_in_pool_cnt = len(database) * CFG_DATABASE_PART / CFG_POOLS_CNT
chain_number = 0
for big_key in database:
	database_pools[chain_number%CFG_POOLS_CNT][big_key] = database[big_key]
	if chain_number > len(database) * CFG_DATABASE_PART:
		break
	chain_number+=1
# print(str(CFG_POOLS_CNT) + ' pools made from ' + str(CFG_DATABASE_PART*100)\
	#   + '% of the database.\nEach pool consists of ' + str(len(database_pools[0])) + " chains")

# output = open("output.txt",)

# for small_key in small_map:
small = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
# merged_set = []
print('=================================\nsmall = ' + small)

def lol(database):
	name_proc = mp.current_process().name
	print(name_proc, len(database))
	return [name_proc]

def getSetForPool(database):
	# 	database = tqdm(database)
	pset = []
	index = 0
	databasesize = len(database)
	for big_key in database:
		# print(index)
		index += 1
		if  index%10 == 0:
			print(mp.current_process().name + ': ' + str(index) + '/' + str(databasesize) + " " + str(len(pset)) + ' found')
		big = database[big_key]
		def calcDP(begin):
			small_len = len(small)
			up_len = math.floor(small_len*CFG_TARGET_SIMILARITY)
			if len(big) < begin + up_len:
				return None, None
			up = big[begin:begin + up_len]
			dp = [[0 for x in range(up_len+1)] for y in range(small_len+1)] # create mat (len(small)+1) X (start_up_len+1)
			# dp = np.zeros((small_len+1, up_len+1))

			for j in range(up_len+1):
				dp[0][j] = j
			for i in range(small_len+1):
				dp[i][0] = i
			for i in range(1, small_len + 1):
				for j in range(1, up_len + 1):
					if small[i-1] == up[j-1]:
						dp[i][j] = dp[i-1][j-1]
					else:
						dp[i][j] = min(dp[i-1][j], dp[i][j-1], dp[i-1][j-1]) + 1

			def addColumn(up_len, begin):
				up_len += 1
				if len(big) < begin + up_len:
					return None

				for i in range(small_len+1):
					if i == 0:
						dp[i].append(up_len)
					elif small[i-1] == big[begin + up_len - 1]:
						dp[i].append(dp[i-1][up_len-1])
					else:
						dp[i].append(min(dp[i-1][up_len], dp[i][up_len-1], dp[i-1][up_len-1]) + 1)

				return up_len

			while dp[small_len][up_len] <= dp[small_len][up_len - 1]:
				new_up_len = addColumn(up_len, begin)
				if new_up_len is not None and dp[small_len][new_up_len] <= dp[small_len][new_up_len - 1]:
					up_len = new_up_len
				else:
					break
			# print(big[begin:begin + up_len] + ' with lev = ' + str(dp[small_len][up_len]))

			return big[begin:begin + up_len], dp[small_len][up_len]

		begin = 0
		while True:
			sub_big, lev = calcDP(begin)
			if lev is None:
				break
			elif lev <= math.floor(len(small)*(1.0 - CFG_TARGET_SIMILARITY)):
				pset.append([lev, sub_big, big_key, begin])
			begin += 1
	return pset
		

if __name__ == '__main__':
	mp.freeze_support()
	with mp.Pool(CFG_POOLS_CNT) as pool:
		list_of_sets = pool.map(getSetForPool, database_pools)
		merged_set = []
		for set1 in list_of_sets:
			merged_set.extend(set1)
		merged_set.sort()
		output = open("output.txt", "w")
		for entry in merged_set:
			output.write(str(entry) + "\n")
		output.close()
		# print(merged_set)