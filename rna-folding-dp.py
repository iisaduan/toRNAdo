def t(rna):
  dp = [[0 for _ in rna] for _ in rna]
  history = [[() for _ in rna] for _ in rna]
  test = []
  for i in range(len(rna)-2, -1, -1):
    for j in range(i+1, len(rna)):
      # LastStep keeps track of max numFolds and possible steps
      max_folds = 0
      possible_folds = []
      # evaluates possble use cases. finds folds for dp algorithm. updates "use" along the way
      for k in range(i+1,j+1):
        compareNum = 0 
        nextFolds = []
        # if you can use it
        if rna[i] == match[rna[k]]:
          nextFolds = [(i,k)]
          # try:
          if i+1 == j:
            compareNum = 1
          elif i+1==k:
            compareNum = 1+dp[k+1][j]
          elif k==j:
            # use = max((use, 1+dp[i+1][j-1]))
            compareNum = 1+dp[i+1][j-1]
          else:
            # use = max((use, 1+dp[i+1][k-1]+dp[k+1][j]))
            compareNum = 1+dp[i+1][k-1]+dp[k+1][j]
          # except:
          #   return (i,j)

          # if the current number of folds is not the maximum
          if max_folds < compareNum: # strictly better
            max_folds = compareNum
            possible_folds = nextFolds
          elif max_folds == compareNum:
            possible_folds.extend(nextFolds)
          else:
            pass # do nothing; the previously found solutions were better
        test.append((i,j,k, possible_folds))
      # reports max of "use" or "lose"
      lose_folds = dp[i+1][j] 
      # try:
      if max_folds < lose_folds: # strictly better
        max_folds = lose_folds
        possible_folds = ["lose " + str(i)]
      elif max_folds == lose_folds:
        possible_folds.append("lose " + str(i))
      else:
        pass # do nothing; the previously found solutions were better
    
      # except: 
        # return lose_folds 

      history[i][j] = possible_folds
      dp[i][j] = max_folds
  return dp, history, test


rna = "AATCAACTA"
dp, history, test = t(rna)
