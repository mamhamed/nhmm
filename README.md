nhmm
====


A "Non-homogeneous Hidden Markov Model" ("NHMM") is a statistical model that extends the Hidden Markov Model (HMM). In a basic hidden Markov model (HMM), there exists of two main processes: 1) an observed or measured process that can be directly measured and 2) a hidden process that cannot been observed directly but it affects the observed process. The observed process (such as user click process in a website) is assumed to be conditionally temporally independent given the hidden process (such as ). The hidden process is assumed to evolve according to a first order Markov chain. In a HMM process state transition probability is fixed with respect to time. A NHMM extends this model by allowing the transition matrix of the hidden states to be dependent on a set of observed covariates . In an example of users click process in a shopping website, clearly the pattern of the click changes by the time of the day. In many applications, the covariates models the underlying seasonality that exists in the hidden process process.

An important application of NHMM is forecasting. Here, the weather states is the hidden states that have a Markov property. In this example, the transition probability between the states depends on the quarter of the year and it almost repeats every year. Hence, the transition probabilities have seasonality factor that is modeled in the NHMM covariates.


Main Contributor:
1- Babak Javid
2- Hamed Firooz
