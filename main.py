from texttable import Texttable
import numpy as np
import math
import re

CFG_TARGET_SIMILARITY = 0.5 # [0.5 .. 1.0]

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

small_map = readChainsMap("lab0.txt")
database = readChainsMap("uniprot_sprot.fasta")

for small_key in small_map:
	small = small_map[small_key]
	set = []
	print('=================================\nsmall = ' + small)
	for big_key in database:
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
				set.append([lev, sub_big, big_key, begin])
			begin += 1
	set.sort()
	print(set)