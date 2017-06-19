#Not completed and not being used yet
#http://crypto.stackexchange.com/questions/25144/super-simple-encryption-of-short-strings
#pseudorandom permutation of 2- or 3-letter initials -> 3-letter pseudonymic initials
#The space is 26*26*27. There's 27 for the last initial because it has 27 possibilities including non-existence
#So take all unique possibilities and randomly shuffle them.

#A matrix 26*26*27. Each entry has a three-letter string. No three-letter string is ever repeated.
#Represent letters and blank by 0 to 26. Or, just have it be a list of the integers from 0 to 26*26*27.
#Then, randomly permute and use modulus to unpack for each string.

allPossibleInitials<- seq(1,26*26*27)
permuted<- sample(allPossibleInitials)
#Either set.seed or gregmisc or 
#http://en.wikipedia.org/wiki/Hash_function#Trivial_hash_function


#Now just use as a look-up table with the string you want to encrypt
initials<- "AOH"
initials<- "ZZ "
numericInitials<- c(0,0,0)
for (i in 1:3) {
	numericInitials[i] <- utf8ToInt(substr(initials,i,i))-64
	if (substr(initials,i,i)==' ')
		numericInitials[i]<- 27
}
idx<-(numericInitials[1]-1)*(26*27) + (numericInitials[2]-1)*27 + numericInitials[3]-1
#idx goes from 0 to 26*26*27-1
encryptedNumeric<- permuted[idx+1]
pseudonymNums<- c("*","*","*")
pseudonymNums[1]<- encryptedNumeric[1] %/% (26*27)
pseudonymNums[2]<- (encryptedNumeric[1] %% (26*27)) %/% (27)
pseudonymNums[3]<- (encryptedNumeric[1] %% (26*27)) %% (27)

letterA = utf8ToInt("A")
pseudonym<- c("*","*","*")
for (i in 1:3) {
	pseudonym[i]<- intToUtf8( letterA + as.numeric(pseudonymNums[i]) )
}
cat("The pseudonymous intials are ",pseudonym)

#decrypt


9*10^2 + 9*10^1 + 9
#hundreds= 567 %/% 100
#tens= (567 %% 100) %/% 10
#ones = (567 %% 100) %% 10

#What is relationship to assigning each subject a number?, S1, S2 etc? 
#Should we avoid 3-letter words?

# Time stamps can easily be manipulated. To avoid detect data corruption and manipulation, I keep lists of MD5 hashes of my crucial data sets. Even if one bit of information is changed, the MD5 hash doesn't match any more. As a precaution, researchers may consider uploading the data declaring its valid MD5 hash somewhere in the paper or supplementary materials. There is excellent freeware for comparing files based on MD5 hash, e.g. MD5check. - Robin Kok
# In a later issue of Mind the Brain, I will describe an horrendous incident in which authors modified a data set available on the web without notice after PPPRs used the data for re-analyses. The authors then used the modified data to redo analyses in an attempt to humiliate the PPPRs with claims that they did not know what they were doing. Fortunately, the PPPRs retained time stamped copies of both data sets. You may like to think that such precautions are unnecessary, but just imagine what critics of PPPR would be saying if they had not saved this evidence.
