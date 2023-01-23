"""
Name:         ensemble
Version:      150212
Author:       Chuan Wang, Lee Chae
Description:  The ensemble module implements ensemble integration algorithms needed by
              E2P2 (Ensemble Enzyme Prediction Pipeline) to provide a final classification
              for a given query protein.

"""

from operator import itemgetter

# Classes
class FinalPredictions():
    def __init__(self, n):
        self.name = n
    def predictions(self):
        self.predictions = []
    def classifiers(self):
        self.classifiers = {}
    def weight(self):
        self.weight = {}

# Functions
def perform_plurality(sequence_id, classifiers, threshold):
    """
    This function reads in a set of sequence objects and a set of classifier objects. For each sequence object,
    it will produce a final enzyme function class prediction via a plurality-rules voting scheme, in which the
    chosen prediction is the one that appears most in the pool of votes. Returns 
    """
    # Instantiate a final prediction object.
    fpred = FinalPredictions(sequence_id)
    fpred.predictions = []
    fpred.classifiers = {}
    final_prediction = ""
    
    # Get all predicted classes for that sequence from all classifiers.
    votes = {}
    for cname in classifiers:
        c = classifiers[cname]
        if sequence_id in c.predictions:
            # Record classifier's predictions for that ID.
            fpred.classifiers[cname] = c.predictions[sequence_id]

            # Tally the predictions.            
            preds = c.predictions[sequence_id]
            for p in preds:
                if votes.has_key(p):
                    votes[p] += 1
                else:
                    votes[p] = 1
                
    # Identify class with the most votes.
    x = len(votes)
    if x == 0:
        # For situations where no classifier has predicted a class for this sequence.
        fpred.predictions.append("NA")
    elif x == 1:
        # For situations where only one vote has been cast.
        for p in votes:
            fpred.predictions.append(p)
    else:
        # Find plurality winner by ranking the classes in order of their vote count.
        sorted_vote = sorted(votes.iteritems(), key=lambda (k,v): (v, k), reverse=True)

        # Get top voted class.
        top_efs = []
        high = sorted_vote[0][1]
        for i in range(len(sorted_vote)):
            vote_count = sorted_vote[i][1]
            if vote_count == high:
                predicted_ef = sorted_vote[i][0]
                top_efs.append(predicted_ef)                    

        # Check for a tie.
        check = len(top_efs)
        if check == 1:
            fpred.predictions.append(sorted_vote[0][0])
        else:
            # No tie-breaker.
            fpred.predictions.append("NA")

    return(fpred)


def perform_majority(sequence_id, classifiers, threshold):
    """
    This function reads in a sequence ID and a set of classifier objects. For each sequence ID,
    it will produce a final enzyme function class prediction via a majority-rules voting scheme, 
    in which the chosen prediction is the one that has the majority (greater than 50%) of votes 
    in the pool.
    """
    # Instantiate a final prediction object.
    fpred = FinalPredictions(sequence_id)
    fpred.predictions = []
    fpred.classifiers = {}

    # Get all predicted classes for that sequence from all classifiers.
    votes = {}
    total_votes = 0
    for cname in classifiers:
        c = classifiers[cname]
        if sequence_id in c.predictions:
            # Record classifier's predictions for that ID.
            fpred.classifiers[cname] = c.predictions[sequence_id]            
            # Tally the predictions.
            preds = c.predictions[sequence_id]
            for p in preds:
                total_votes += 1
                if votes.has_key(p):
                    votes[p] += 1
                else:
                    votes[p] = 1
    
    # Identify class with the most votes.
    x = len(votes)
    if x == 0:
        # For situations where no classifier has predicted a class for this sequence.
        fpred.predictions.append("NA")
    elif x == 1:
        # For situations where only one vote has been cast.
        for p in votes:
            fpred.predictions.append(p)
    else:
        # Find plurality winner by ranking the classes in order of their vote count.
        sorted_vote = sorted(votes.iteritems(), key=lambda (k,v): (v, k), reverse=True)

        # Identify top voted class.
        top_efs = []
        high = sorted_vote[0][1]
        
        # Check to make sure there is a majority winner.
        majority_needed = (total_votes/2)
        if high > majority_needed:
            for i in range(len(sorted_vote)):
                vote_count = sorted_vote[i][1]
                if vote_count == high:
                    predicted_ef = sorted_vote[i][0]
                    top_efs.append(predicted_ef)

            # Check for a tie.
            check = len(top_efs)
            if check == 1:
                fpred.predictions.append(sorted_vote[0][0])
            else:
                # No tie-breaker.
                fpred.predictions.append("NA")
        else:
            fpred.predictions.append("NA")
            
    return(fpred)
        
def perform_max_weight(sequence_id, classifiers, threshold):
    """
    This function reads in a sequence ID and a set of classifier objects. For each sequence ID,    
    it will produce a final enzyme function class prediction via a maximum weight voting scheme, 
    in which the chosen prediction is the one that has the greatest weight associated with it. 
    Weights are performance measures obtained from training data. Note that a numerical threshold 
    can be invoked that would allow other votes within a certain range of the winning weight to also 
    be included in the final prediction. This threshold allows for the presence of multi-function
    enzymes.
    """
    # Instantiate a final prediction object.
    fpred = FinalPredictions(sequence_id)
    fpred.predictions = []
    fpred.classifiers = {}
        
    # Get all predicted classes for that sequence from all classifiers.
    votes = []
    total_votes = 0
    for cname in classifiers:
        c = classifiers[cname]
        if sequence_id in c.predictions:
            temp = []
            # Record classifier's predictions and weights for that ID for voting and for output.
            for ef_class in c.predictions[sequence_id]:
                try:
                    ef_weight = c.weights[ef_class]
                except:
                    ef_weight = "0.000"
                # Record vote as tuple, as this will allow one EF class to have more than one weight, depending
                # on if more than one classifier called it.
                vote = (cname, ef_class, float(ef_weight))
                votes.append(vote)
                # The classifiers attribute will be used in outputting the full results. 
                entry = "%s (%s)" % (ef_class, ef_weight)
                temp.append(entry)
            fpred.classifiers[cname] = temp
 
    # Check to see if any votes were recorded.
    x = len(votes)
    if x < 1:
        fpred.predictions.append("NA")
    else:
        # Sort the votes to find the highest weighted vote.
        sorted_weights = sorted(votes, key=itemgetter(2), reverse=True)
        high_weight = sorted_weights[0][2]
        high_class = sorted_weights[0][1]
        entry = "%s (%s)" % (high_class, str(high_weight))
        fpred.predictions.append(entry)

        # Condense the votes so that each EF class appears only once, with its highest weight.
        condensed_votes = {}
        for vote in votes:
            ef_class = vote[1]
            weight = float(vote[2])
            if ef_class in condensed_votes:
                if weight > condensed_votes[ef_class]:
                    condensed_votes[ef_class] = weight
            else:
                condensed_votes[ef_class] = weight

        # Iterate through the votes to find all those that pass the threshold.
        t = high_weight - ( float(threshold)/float(100) )
        for ef_class in condensed_votes:
            if ef_class != high_class:
                weight = condensed_votes[ef_class]
                if float(weight) >= t:
                    entry = "%s (%s)" % (ef_class, weight)
                    fpred.predictions.append(entry)
    return(fpred)

def perform_max_weight_percent_threshold(sequence_id, classifiers, threshold):
    """
    This function reads in a sequence ID and a set of classifier objects. For each sequence ID,    
    it will produce a final enzyme function class prediction via a maximum weight voting scheme, 
    in which the chosen prediction is the one that has the greatest weight associated with it. 
    Weights are performance measures obtained from training data. Note that a numerical threshold 
    can be invoked that would allow other votes within a certain range of the winning weight to also 
    be included in the final prediction. This threshold allows for the presence of multi-function
    enzymes.
    """
    # Instantiate a final prediction object.
    fpred = FinalPredictions(sequence_id)
    fpred.predictions = []
    fpred.classifiers = {}
        
    # Get all predicted classes for that sequence from all classifiers.
    votes = []
    total_votes = 0
    for cname in classifiers:
        c = classifiers[cname]
        if sequence_id in c.predictions:
            temp = []
            # Record classifier's predictions and weights for that ID for voting and for output.
            for ef_class in c.predictions[sequence_id]:
                try:
                    ef_weight = c.weights[ef_class]
                except:
                    ef_weight = "0.000"
                # Record vote as tuple, as this will allow one EF class to have more than one weight, depending
                # on if more than one classifier called it.
                vote = (cname, ef_class, float(ef_weight))
                votes.append(vote)
                # The classifiers attribute will be used in outputting the full results. 
                entry = "%s (%s)" % (ef_class, ef_weight)
                temp.append(entry)
            fpred.classifiers[cname] = temp
 
    # Check to see if any votes were recorded.
    x = len(votes)
    if x < 1:
        fpred.predictions.append("NA")
    else:
        # Sort the votes to find the highest weighted vote.
        sorted_weights = sorted(votes, key=itemgetter(2), reverse=True)
        high_weight = sorted_weights[0][2]
        high_class = sorted_weights[0][1]
        entry = "%s (%s)" % (high_class, str(high_weight))
        fpred.predictions.append(entry)

        # Condense the votes so that each EF class appears only once, with its highest weight.
        condensed_votes = {}
        for vote in votes:
            ef_class = vote[1]
            weight = float(vote[2])
            if ef_class in condensed_votes:
                if weight > condensed_votes[ef_class]:
                    condensed_votes[ef_class] = weight
            else:
                condensed_votes[ef_class] = weight

        # Iterate through the votes to find all those that pass the threshold.
        t = high_weight * ( float(1 - (float(threshold)/float(100)) ) )
        for ef_class in condensed_votes:
            if ef_class != high_class:
                weight = condensed_votes[ef_class]
                if float(weight) >= t:
                    entry = "%s (%s)" % (ef_class, weight)
                    fpred.predictions.append(entry)
    return(fpred)
        
def perform_max_weight_absolute_threshold(sequence_id, classifiers, threshold):
    """
    This function reads in a sequence ID and a set of classifier objects. For each sequence ID,    
    it will produce a final enzyme function class prediction via a maximum weight voting scheme, 
    in which the chosen prediction is the one that has the greatest weight associated with it. 
    Weights are performance measures obtained from training data. Note that a numerical threshold 
    can be invoked that would allow other votes within a certain range of the winning weight to also 
    be included in the final prediction. This threshold allows for the presence of multi-function
    enzymes.
    """
    # Instantiate a final prediction object.
    fpred = FinalPredictions(sequence_id)
    fpred.predictions = []
    fpred.classifiers = {}
        
    # Get all predicted classes for that sequence from all classifiers.
    votes = []
    total_votes = 0
    for cname in classifiers:
        c = classifiers[cname]
        if sequence_id in c.predictions:
            temp = []
            # Record classifier's predictions and weights for that ID for voting and for output.
            for ef_class in c.predictions[sequence_id]:
                try:
                    ef_weight = c.weights[ef_class]
                except:
                    ef_weight = "0.000"
                # Record vote as tuple, as this will allow one EF class to have more than one weight, depending
                # on if more than one classifier called it.
                vote = (cname, ef_class, float(ef_weight))
                votes.append(vote)
                # The classifiers attribute will be used in outputting the full results. 
                entry = "%s (%s)" % (ef_class, ef_weight)
                temp.append(entry)
            fpred.classifiers[cname] = temp
 
    # Check to see if any votes were recorded.
    x = len(votes)
    if x < 1:
        fpred.predictions.append("NA")
    else:
        # Sort the votes to find the highest weighted vote.
        sorted_weights = sorted(votes, key=itemgetter(2), reverse=True)
        high_weight = sorted_weights[0][2]
        high_class = sorted_weights[0][1]
        entry = "%s (%s)" % (high_class, str(high_weight))
        fpred.predictions.append(entry)

        # Condense the votes so that each EF class appears only once, with its highest weight.
        condensed_votes = {}
        for vote in votes:
            ef_class = vote[1]
            weight = float(vote[2])
            if ef_class in condensed_votes:
                if weight > condensed_votes[ef_class]:
                    condensed_votes[ef_class] = weight
            else:
                condensed_votes[ef_class] = weight

        # Iterate through the votes to find all those that pass the threshold.
        t = float(high_weight - threshold)
        # Check to make sure the threshold is not negative.
        if t < 0.0:
            t = 0.0
        for ef_class in condensed_votes:
            if ef_class != high_class:
                weight = condensed_votes[ef_class]
                if float(weight) >= t:
                    entry = "%s (%s)" % (ef_class, weight)
                    fpred.predictions.append(entry)
    return(fpred)
        
def perform_avg_weight(sequence_id, classifiers, threshold):
    """
    This function reads in a sequence ID and a set of classifier objects. For each sequence ID,
    it will produce a final enzyme function class prediction via an average weight voting scheme, 
    in which the chosen prediction is the one that has the greatest average weight associated with 
    it. Weights are performance measures obtained from training data. Note that a numerical threshold 
    can be invoked that would allow other votes within a certain range of the winning weight to also 
    be included in the final prediction. This threshold allows for the presence of multi-function
    enzymes.
    """
    # Instantiate a final prediction object.
    fpred = FinalPredictions(sequence_id)
    fpred.predictions = []
    fpred.classifiers = {}
        
    # Get all predicted classes for that sequence from all classifiers.
    votes = {}
    total_votes = 0
    for cname in classifiers:
        c = classifiers[cname]
        if sequence_id in c.predictions:
            temp = []
            # Record classifier's predictions and weights for that ID for voting and for output.
            for ef_class in c.predictions[sequence_id]:
                try:
                    ef_weight = c.weights[ef_class]
                except:
                    ef_weight = "0.000"
                # For each vote, add its weight. Keep track of all votes cast.
                total_votes += 1
                if votes.has_key(ef_class):
                    votes[ef_class] += float(ef_weight)
                else:
                    votes[ef_class] = float(ef_weight)
                # The classifiers attribute will be used in outputting the full results. 
                entry = "%s (%s)" % (ef_class, ef_weight)
                temp.append(entry)
            fpred.classifiers[cname] = temp


    # Add routines to find no-vote and single-vote events.
    # Check to see if any votes were recorded.
    x = len(votes)
    if x < 1:
        fpred.predictions.append("NA")
    else:
        # For each vote, find its average weight.
        avg_weights = {}
        for ef_class in votes:
            total_weight = votes[ef_class]
            avg_weights[ef_class] = float(total_weight/total_votes)
         
        # Sort the votes from highest average weight to lowest.
        sorted_avg_weights = sorted(avg_weights.iteritems(), key=lambda (k,v): (v, k), reverse=True)   
        high_class = sorted_avg_weights[0][0]
        high_weight = sorted_avg_weights[0][1]
        entry = "%s (%s)" % (high_class, str(high_weight))
        fpred.predictions.append(entry)
     
        # Iterate through the votes to find all votes with average weights within the threshold
        # of the top weight.
        t = high_weight - ( float(threshold)/float(100) )
        for class_weight_pair in sorted_avg_weights[1:]: #The top class has already been recorded.
            ef_class = class_weight_pair[0]
            weight = class_weight_pair[1]
            if float(weight) >= t:
                entry = "%s (%s)" % (ef_class, weight)
                fpred.predictions.append(entry)
    return(fpred)

def perform_avg_weight_percent_threshold(sequence_id, classifiers, threshold):
    """
    This function reads in a sequence ID and a set of classifier objects. For each sequence ID,
    it will produce a final enzyme function class prediction via an average weight voting scheme, 
    in which the chosen prediction is the one that has the greatest average weight associated with 
    it. Weights are performance measures obtained from training data. Note that a numerical threshold 
    can be invoked that would allow other votes within a certain range of the winning weight to also 
    be included in the final prediction. This threshold allows for the presence of multi-function
    enzymes.
    """
    # Instantiate a final prediction object.
    fpred = FinalPredictions(sequence_id)
    fpred.predictions = []
    fpred.classifiers = {}
        
    # Get all predicted classes for that sequence from all classifiers.
    votes = {}
    total_votes = 0
    for cname in classifiers:
        c = classifiers[cname]
        if sequence_id in c.predictions:
            temp = []
            # Record classifier's predictions and weights for that ID for voting and for output.
            for ef_class in c.predictions[sequence_id]:
                try:
                    ef_weight = c.weights[ef_class]
                except:
                    ef_weight = "0.000"
                # For each vote, add its weight. Keep track of all votes cast.
                total_votes += 1
                if votes.has_key(ef_class):
                    votes[ef_class] += float(ef_weight)
                else:
                    votes[ef_class] = float(ef_weight)
                # The classifiers attribute will be used in outputting the full results. 
                entry = "%s (%s)" % (ef_class, ef_weight)
                temp.append(entry)
            fpred.classifiers[cname] = temp


    # Add routines to find no-vote and single-vote events.
    # Check to see if any votes were recorded.
    x = len(votes)
    if x < 1:
        fpred.predictions.append("NA")
    else:
        # For each vote, find its average weight.
        avg_weights = {}
        for ef_class in votes:
            total_weight = votes[ef_class]
            avg_weights[ef_class] = float(total_weight/total_votes)
         
        # Sort the votes from highest average weight to lowest.
        sorted_avg_weights = sorted(avg_weights.iteritems(), key=lambda (k,v): (v, k), reverse=True)   
        high_class = sorted_avg_weights[0][0]
        high_weight = sorted_avg_weights[0][1]
        entry = "%s (%s)" % (high_class, str(high_weight))
        fpred.predictions.append(entry)
     
        # Iterate through the votes to find all votes with average weights within the threshold
        # of the top weight.
        t = high_weight * ( float(1 - (float(threshold)/float(100)) ) ) 
        for class_weight_pair in sorted_avg_weights[1:]: #The top class has already been recorded.
            ef_class = class_weight_pair[0]
            weight = class_weight_pair[1]
            if float(weight) >= t:
                entry = "%s (%s)" % (ef_class, weight)
                fpred.predictions.append(entry)
    return(fpred)

def perform_avg_weight_absolute_threshold(sequence_id, classifiers, threshold):
    """
    This function reads in a sequence ID and a set of classifier objects. For each sequence ID,
    it will produce a final enzyme function class prediction via an average weight voting scheme, 
    in which the chosen prediction is the one that has the greatest average weight associated with 
    it. Weights are performance measures obtained from training data. Note that a numerical threshold 
    can be invoked that would allow other votes within a certain range of the winning weight to also 
    be included in the final prediction. This threshold allows for the presence of multi-function
    enzymes.
    """
    # Instantiate a final prediction object.
    fpred = FinalPredictions(sequence_id)
    fpred.predictions = []
    fpred.classifiers = {}
        
    # Get all predicted classes for that sequence from all classifiers.
    votes = {}
    total_votes = 0
    for cname in classifiers:
        c = classifiers[cname]
        if sequence_id in c.predictions:
            temp = []
            # Record classifier's predictions and weights for that ID for voting and for output.
            for ef_class in c.predictions[sequence_id]:
                try:
                    ef_weight = c.weights[ef_class]
                except:
                    ef_weight = "0.000"
                # For each vote, add its weight. Keep track of all votes cast.
                total_votes += 1
                if votes.has_key(ef_class):
                    votes[ef_class] += float(ef_weight)
                else:
                    votes[ef_class] = float(ef_weight)
                # The classifiers attribute will be used in outputting the full results. 
                entry = "%s (%s)" % (ef_class, ef_weight)
                temp.append(entry)
            fpred.classifiers[cname] = temp


    # Add routines to find no-vote and single-vote events.
    # Check to see if any votes were recorded.
    x = len(votes)
    if x < 1:
        fpred.predictions.append("NA")
    else:
        # For each vote, find its average weight.
        avg_weights = {}
        for ef_class in votes:
            total_weight = votes[ef_class]
            avg_weights[ef_class] = float(total_weight/total_votes)
         
        # Sort the votes from highest average weight to lowest.
        sorted_avg_weights = sorted(avg_weights.iteritems(), key=lambda (k,v): (v, k), reverse=True)   
        high_class = sorted_avg_weights[0][0]
        high_weight = sorted_avg_weights[0][1]
        entry = "%s (%s)" % (high_class, str(high_weight))
        fpred.predictions.append(entry)
     
        # Iterate through the votes to find all votes with average weights within the threshold
        # of the top weight.
        t = float(high_weight - threshold)
        if t < 0.0:
            t = 0.0
        for class_weight_pair in sorted_avg_weights[1:]: #The top class has already been recorded.
            ef_class = class_weight_pair[0]
            weight = class_weight_pair[1]
            if float(weight) >= t:
                entry = "%s (%s)" % (ef_class, weight)
                fpred.predictions.append(entry)
    return(fpred)

def perform_fixed_avg_weight_absolute_threshold(sequence_id, classifiers, threshold):
    """
    This function reads in a sequence ID and a set of classifier objects. For each sequence ID,
    it will produce a final enzyme function class prediction via an average weight voting scheme, 
    in which the chosen prediction is the one that has the greatest average weight associated with 
    it. Weights are performance measures obtained from training data. Note that a numerical threshold 
    can be invoked that would allow other votes within a certain range of the winning weight to also 
    be included in the final prediction. This threshold allows for the presence of multi-function
    enzymes.
    """
    # Instantiate a final prediction object.
    fpred = FinalPredictions(sequence_id)
    fpred.predictions = []
    fpred.classifiers = {}
        
    # Get all predicted classes for that sequence from all classifiers.
    votes = {}
    total_votes = {}
    for cname in classifiers:
        c = classifiers[cname]
        if sequence_id in c.predictions:
            temp = []
            # Record classifier's predictions and weights for that ID for voting and for output.
            for ef_class in c.predictions[sequence_id]:
                try:
                    ef_weight = c.weights[ef_class]
                except:
                    ef_weight = "0.000"
                # For each vote, add its weight. Keep track of all votes cast.
                if total_votes.has_key(ef_class):
                    total_votes[ef_class] += 1
                else:
                    total_votes[ef_class] = 1
                if votes.has_key(ef_class):
                    votes[ef_class] += float(ef_weight)
                else:
                    votes[ef_class] = float(ef_weight)
                # The classifiers attribute will be used in outputting the full results. 
                entry = "%s (%s)" % (ef_class, ef_weight)
                temp.append(entry)
            fpred.classifiers[cname] = temp


    # Add routines to find no-vote and single-vote events.
    # Check to see if any votes were recorded.
    x = len(votes)
    if x < 1:
        fpred.predictions.append("NA")
    else:
        # For each vote, find its average weight.
        avg_weights = {}
        for ef_class in votes:
            total_weight = votes[ef_class]
            avg_weights[ef_class] = float(total_weight/total_votes[ef_class])
         
        # Sort the votes from highest average weight to lowest.
        sorted_avg_weights = sorted(avg_weights.iteritems(), key=lambda (k,v): (v, k), reverse=True)   
        high_class = sorted_avg_weights[0][0]
        high_weight = sorted_avg_weights[0][1]
        entry = "%s (%s)" % (high_class, str(high_weight))
        fpred.predictions.append(entry)
     
        # Iterate through the votes to find all votes with average weights within the threshold
        # of the top weight.
        t = float(high_weight - threshold)
        if t < 0.0:
            t = 0.0
        for class_weight_pair in sorted_avg_weights[1:]: #The top class has already been recorded.
            ef_class = class_weight_pair[0]
            weight = class_weight_pair[1]
            if float(weight) >= t:
                entry = "%s (%s)" % (ef_class, weight)
                fpred.predictions.append(entry)
    return(fpred)

def perform_binary_avg_weight_absolute_threshold(sequence_id, classifiers, threshold):
    """
    This function reads in a sequence ID and a set of classifier objects. For each sequence ID,
    it will produce a final enzyme function class prediction via an binary average weight voting scheme, 
    in which the chosen prediction is the one that has the greatest average weight associated with 
    it. Weights are performance measures obtained from training data. Note that a numerical threshold 
    can be invoked that would allow other votes within a certain range of the winning weight to also 
    be included in the final prediction. This threshold allows for the presence of multi-function
    enzymes.
    """
    # Instantiate a final prediction object.
    fpred = FinalPredictions(sequence_id)
    fpred.predictions = []
    fpred.classifiers = {}
        
    # Get all predicted classes for that sequence from all classifiers.
    votes = {}
    total_votes = 0
    all_ef = {}
    for cname in classifiers:
        c = classifiers[cname]
        if sequence_id in c.predictions:
            temp = []
            # Record classifier's predictions and weights for that ID for voting and for output.
            for ef_class in c.predictions[sequence_id]:
                try:
                    ef_weight = c.weights[ef_class]
                except:
                    ef_weight = "0.000"
                # For each vote, add its weight. Keep track of all votes cast.
                total_votes += 1
                all_ef[ef_class] = 1
                if votes.has_key(ef_class):
                    votes[ef_class] += float(ef_weight)
                else:
                    votes[ef_class] = float(ef_weight)
                # The classifiers attribute will be used in outputting the full results. 
                entry = "%s (%s)" % (ef_class, ef_weight)
                temp.append(entry)
            fpred.classifiers[cname] = temp
    
    for ef in all_ef:
        for cname in classifiers:
            c = classifiers[cname]
            if sequence_id in c.predictions:
                seqpred = {}
                for ef_class in c.predictions[sequence_id]:
                    seqpred[ef_class] = 1
                if not seqpred.has_key(ef):
                    try:
                        ef_weight = c.weights[ef_class]
                    except:
                        ef_weight = "0.000"
                    if votes.has_key(ef_class):
                        votes[ef_class] -= float(ef_weight)
                    else:
                        votes[ef_class] = -float(ef_weight)
            else:
                try:
                    ef_weight = c.weights[ef_class]
                except:
                    ef_weight = "0.000"
                if votes.has_key(ef_class):
                    votes[ef_class] -= float(ef_weight)
                else:
                    votes[ef_class] = 0-float(ef_weight)
    

    # Add routines to find no-vote and single-vote events.
    # Check to see if any votes were recorded.
    x = len(votes)
    if x < 1:
        fpred.predictions.append("NA")
    else:
        # For each vote, find its average weight.
        avg_weights = {}
        for ef_class in votes:
            total_weight = votes[ef_class]
            avg_weights[ef_class] = float(total_weight/2)
         
        # Sort the votes from highest average weight to lowest.
        sorted_avg_weights = sorted(avg_weights.iteritems(), key=lambda (k,v): (v, k), reverse=True)   
        high_class = sorted_avg_weights[0][0]
        high_weight = sorted_avg_weights[0][1]
        entry = "%s (%s)" % (high_class, str(high_weight))
        fpred.predictions.append(entry)
     
        # Iterate through the votes to find all votes with average weights within the threshold
        # of the top weight.
        t = float(high_weight - threshold)
        if t < 0.0:
            t = 0.0
        for class_weight_pair in sorted_avg_weights[1:]: #The top class has already been recorded.
            ef_class = class_weight_pair[0]
            weight = class_weight_pair[1]
            if float(weight) >= t:
                entry = "%s (%s)" % (ef_class, weight)
                fpred.predictions.append(entry)
    return(fpred)

def perform_fixed_avg_weight_percent_threshold(sequence_id, classifiers, threshold):
    """
    This function reads in a sequence ID and a set of classifier objects. For each sequence ID,
    it will produce a final enzyme function class prediction via an average weight voting scheme, 
    in which the chosen prediction is the one that has the greatest average weight associated with 
    it. Weights are performance measures obtained from training data. Note that a numerical threshold 
    can be invoked that would allow other votes within a certain range of the winning weight to also 
    be included in the final prediction. This threshold allows for the presence of multi-function
    enzymes.
    """
    # Instantiate a final prediction object.
    fpred = FinalPredictions(sequence_id)
    fpred.predictions = []
    fpred.classifiers = {}
        
    # Get all predicted classes for that sequence from all classifiers.
    votes = {}
    total_votes = {}
    for cname in classifiers:
        c = classifiers[cname]
        if sequence_id in c.predictions:
            temp = []
            # Record classifier's predictions and weights for that ID for voting and for output.
            for ef_class in c.predictions[sequence_id]:
                try:
                    ef_weight = c.weights[ef_class]
                except:
                    ef_weight = "0.000"
                # For each vote, add its weight. Keep track of all votes cast.
                if total_votes.has_key(ef_class):
                    total_votes[ef_class] += 1
                else:
                    total_votes[ef_class] = 1
                if votes.has_key(ef_class):
                    votes[ef_class] += float(ef_weight)
                else:
                    votes[ef_class] = float(ef_weight)
                # The classifiers attribute will be used in outputting the full results. 
                entry = "%s (%s)" % (ef_class, ef_weight)
                temp.append(entry)
            fpred.classifiers[cname] = temp


    # Add routines to find no-vote and single-vote events.
    # Check to see if any votes were recorded.
    x = len(votes)
    if x < 1:
        fpred.predictions.append("NA")
    else:
        # For each vote, find its average weight.
        avg_weights = {}
        for ef_class in votes:
            total_weight = votes[ef_class]
            avg_weights[ef_class] = float(total_weight/total_votes[ef_class])
         
        # Sort the votes from highest average weight to lowest.
        sorted_avg_weights = sorted(avg_weights.iteritems(), key=lambda (k,v): (v, k), reverse=True)   
        high_class = sorted_avg_weights[0][0]
        high_weight = sorted_avg_weights[0][1]
        entry = "%s (%s)" % (high_class, str(high_weight))
        fpred.predictions.append(entry)
     
        # Iterate through the votes to find all votes with average weights within the threshold
        # of the top weight.
        t = high_weight * ( float(1 - (float(threshold)/float(100)) ) ) 
        if t < 0.0:
            t = 0.0
        for class_weight_pair in sorted_avg_weights[1:]: #The top class has already been recorded.
            ef_class = class_weight_pair[0]
            weight = class_weight_pair[1]
            if float(weight) >= t:
                entry = "%s (%s)" % (ef_class, weight)
                fpred.predictions.append(entry)
    return(fpred)

def perform_binary_avg_weight_percent_threshold(sequence_id, classifiers, threshold):
    """
    This function reads in a sequence ID and a set of classifier objects. For each sequence ID,
    it will produce a final enzyme function class prediction via an binary average weight voting scheme, 
    in which the chosen prediction is the one that has the greatest average weight associated with 
    it. Weights are performance measures obtained from training data. Note that a numerical threshold 
    can be invoked that would allow other votes within a certain range of the winning weight to also 
    be included in the final prediction. This threshold allows for the presence of multi-function
    enzymes.
    """
    # Instantiate a final prediction object.
    fpred = FinalPredictions(sequence_id)
    fpred.predictions = []
    fpred.classifiers = {}
        
    # Get all predicted classes for that sequence from all classifiers.
    votes = {}
    total_votes = 0
    all_ef = {}
    for cname in classifiers:
        c = classifiers[cname]
        if sequence_id in c.predictions:
            temp = []
            # Record classifier's predictions and weights for that ID for voting and for output.
            for ef_class in c.predictions[sequence_id]:
                try:
                    ef_weight = c.weights[ef_class]
                except:
                    ef_weight = "0.000"
                # For each vote, add its weight. Keep track of all votes cast.
                total_votes += 1
                all_ef[ef_class] = 1
                if votes.has_key(ef_class):
                    votes[ef_class] += float(ef_weight)
                else:
                    votes[ef_class] = float(ef_weight)
                # The classifiers attribute will be used in outputting the full results. 
                entry = "%s (%s)" % (ef_class, ef_weight)
                temp.append(entry)
            fpred.classifiers[cname] = temp
    
    for ef in all_ef:
        for cname in classifiers:
            c = classifiers[cname]
            if sequence_id in c.predictions:
                seqpred = {}
                for ef_class in c.predictions[sequence_id]:
                    seqpred[ef_class] = 1
                if not seqpred.has_key(ef):
                    try:
                        ef_weight = c.weights[ef_class]
                    except:
                        ef_weight = "0.000"
                    if votes.has_key(ef_class):
                        votes[ef_class] -= float(ef_weight)
                    else:
                        votes[ef_class] = -float(ef_weight)
            else:
                try:
                    ef_weight = c.weights[ef_class]
                except:
                    ef_weight = "0.000"
                if votes.has_key(ef_class):
                    votes[ef_class] -= float(ef_weight)
                else:
                    votes[ef_class] = -float(ef_weight)
    

    # Add routines to find no-vote and single-vote events.
    # Check to see if any votes were recorded.
    x = len(votes)
    if x < 1:
        fpred.predictions.append("NA")
    else:
        # For each vote, find its average weight.
        avg_weights = {}
        for ef_class in votes:
            total_weight = votes[ef_class]
            avg_weights[ef_class] = float(total_weight/2)
         
        # Sort the votes from highest average weight to lowest.
        sorted_avg_weights = sorted(avg_weights.iteritems(), key=lambda (k,v): (v, k), reverse=True)   
        high_class = sorted_avg_weights[0][0]
        high_weight = sorted_avg_weights[0][1]
        entry = "%s (%s)" % (high_class, str(high_weight))
        fpred.predictions.append(entry)
     
        # Iterate through the votes to find all votes with average weights within the threshold
        # of the top weight.
        t = high_weight * ( float(1 - (float(threshold)/float(100)) ) ) 
        if t < 0.0:
            t = 0.0
        for class_weight_pair in sorted_avg_weights[1:]: #The top class has already been recorded.
            ef_class = class_weight_pair[0]
            weight = class_weight_pair[1]
            if float(weight) >= t:
                entry = "%s (%s)" % (ef_class, weight)
                fpred.predictions.append(entry)
    return(fpred)
