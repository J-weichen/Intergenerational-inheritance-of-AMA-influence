#! python2

window = 1000000 #1M
fout = open("1.1M.bed", "w")

fin = open("/mnt/data/chenwei/qinmen_BR/0.meng_script/1.meth_script/CNV_1Mb_script/hg38.len.pure")
for line in fin:
	if line.startswith("chrM"):
		continue
	col = line.strip().split()
	chrom = col[0]
	chrlen = int(col[1])
	for i in xrange(1, chrlen+1, window):
		fout.write(chrom+"\t"+str(i)+"\t"+str(i+window-1)+"\n")
fin.close()
fout.close()


