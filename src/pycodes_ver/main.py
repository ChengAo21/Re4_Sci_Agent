from langgraph.graph import StateGraph, START, END
from IPython.display import display, Image
from schemas import State
from nodes import llm_call_c1, llm_call_p2, run_codes, llm_call_r3, decide_next_step
from problems import burgers_problem, cavity_flow_problem, sod_shock_problem, longTm_NS_problem, poison_problem, helmholtz_problem, cylinder_problem, hilbert_mat_problem, Keyhole_problem
from utils import save_response_toMD

# Build Graph
agent_builder = StateGraph(State)
agent_builder.add_node("consultant", llm_call_c1)
agent_builder.add_node("programmer", llm_call_p2)
agent_builder.add_node("executor", run_codes)
agent_builder.add_node("reviewer", llm_call_r3)


agent_builder.add_edge(START, "consultant")
agent_builder.add_edge("consultant", "programmer")
agent_builder.add_edge("programmer", "executor")

agent_builder.add_edge("executor", "reviewer")

agent_builder.add_conditional_edges(
    "reviewer",
    decide_next_step,
    {
        "Accepted": END,
        "Revised": "programmer", 
    },
)


agent_graph = agent_builder.compile()

# Visualize Graph if needed
# display(Image(agent_graph.get_graph().draw_mermaid_png()))


# Run Example
if __name__ == "__main__":
    
    # 默认 burgers_problem
    selected_problem = burgers_problem
    
    initial_state = {
        "prob_todo": selected_problem,
    }
    
    final_state = agent_graph.invoke(initial_state)

    save_response_toMD(final_state)