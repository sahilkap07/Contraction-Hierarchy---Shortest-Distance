# Contraction-Hierarchy---Shortest-Distance

## Problem Inroduction :
Dijkstra's and A* algorithms solve the shortest path problem by iteratively relaxing the next closest node with the smallest estimated distance. Contraction hierarchy instead solves the same problem in 2 steps : Preprocessing and Query answering. 
This can reduce the query time to some milliseconds from a few seconds for the query handling step. The idea is to contract the nodes in the reverse order of their importance. Given a graph with nodes and edges connecting them, the problem requires to preprocess the graph and then answer the subsequent queries using the preprocessed graph in an efficient manner.

## Underlying Concept :
The concept is based on preprocessing of the graph to answer the queries in an efficient way. The preprocessing stage includes calculation of importance for each node such that nodes with least importance are contracted first. 
Importance is generally a function of four values calculated in preprocssing stage :
1. Edge Difference (ed) - Difference between the added edges and the already existing edges.
2. Contracted Neighbours (cn) - Number of neighbours contracted.
3. Shortcut Cover (sc) - Number of nodes affected as a result of contraction.
4. Node Level (nl) - Maximum of the node level and the level of the contracted neighbour.

Importance is then - (alpha * ed) + (beta * cn) + (gamma * sc) + (psi * nl), 
                               where alpha, beta, gamma and psi are heuristic multiplier values.
                               
The importance of nodes keep on changing during the extraction process and needs to be recomputed continuously.

Once the nodes are contracted and ranked, queries can then be answered by running a bidirectional dijkstra's from both the source and target nodes till the estimated path length exceeds the already calculated path length till that point of time.

## Input Format : 
The first line will consist of two integers n and m - number of nodes and edges in the graph. Each of the subsequent m lines will contain three integers each, namely u, v, c - initial node number, final node and the weight of that node. The program will then preprocess the graph and print "Ready" when done with preprocessing. Finally, the query part will include inputting an integer t - number of queries and following t lines with two integers each, u and v - initial and final node for the calculation of shortest path.

## Constraints : 
1 ≤ 𝑛 ≤ 110 000; 1 ≤ 𝑚 ≤ 250 000; 1 ≤ 𝑢, 𝑣 ≤ 𝑛; 1 ≤ 𝑙 ≤ 200 000; 1 ≤ 𝑞 ≤ 10 000. It is guaranteed that the correct distances are less than 1 000 000 000.

## Output Format :
The output will have t lines - number of queries asnwering the shortest path length between the nodes.
