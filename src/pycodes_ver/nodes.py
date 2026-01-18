import os
from IPython.display import display, Markdown
from config import MAX_ITERATIONS, MULTI_MODAL_REVIEW, ALTERNATIVE_NUM, MAX_COMMENTS, ConsModule_ON, CUT_OUTEXT, PROG_name, Cons_name, Rev_name
from schemas import State, Consultant_module, Progarmmer_module, Reviewer_module, Reivew_multimodal_image
from utils import collapsible_mardown, execute_code_tool, load_image_tool, encode_image_tool
from llm_core import llm_call_

def llm_call_c1(state: State):
    display(Markdown(f"### 🛌🏻 pdeAgent with Multi-Modal Review"))
    collapsible_mardown("Configs 🏷️", f"MaxCount: {MAX_ITERATIONS}\nImageReview: {MULTI_MODAL_REVIEW}\nPlanNum: {ALTERNATIVE_NUM}\nMaxRevComs: {MAX_COMMENTS}\nConsModule_ON: {ConsModule_ON}\nCUT_OUTEXT: {CUT_OUTEXT}")
    if ConsModule_ON:
        display(Markdown("### 🧠 Consultant:"))
        sys_msg = (
            "You are a distinguished Mathematical Consultant and Numerical Analyst.\n"
            "[Task] Expand the context of a given scientific problem and generate multiple alternative solution options.\n"        
            "[Goal] Accurately dissect the problem by preserving its original meaning, explicitly identifying the underlying mathematical and numerical challenges, and formulating rigorous and well-structured solution strategies."        
        )

        config_dict = {
            "llm_name": Cons_name,
            "system_content": sys_msg,
            "human_content": f"[Problem Statement]:\n{state['prob_todo']}",
        }
        cons_response = llm_call_(config_dict, Consultant_module)

        plans_text = "\n".join(
            f"{i+1}. {pl.solu_name}:\n{pl.content}\n\n" for i,pl in enumerate(cons_response.solution_plans)
        )
        dict_text = {pl.solu_name:pl.content for pl in cons_response.solution_plans}
        # No need to store titles    
        plans_title = "\n".join(dict_text.keys())
        # print(f"Multiple Plans:\n\n{plans_title}")
        collapsible_mardown(f"Multiple Plans:", plans_text)

        return {
            "expanded_prob":cons_response.expanded_prob,
            "solution_plans": cons_response.solution_plans,
            "iteration_count": 0,
        }
    else:
        return {"iteration_count": 0,
                "expanded_prob": "",
                "solution_plans": [],}


def llm_call_p2(state: State):
    display(Markdown("### 👨‍💻 Programmer Generating:"))
    sys_msg_coder = ("You are an expert Python Programmer specialized in scientific computing.\n"
                )
    
    human_msg_coder = (
                f"[Problem Statement]:\n{state['prob_todo']}\n"
                )

    if ConsModule_ON:
        plans_text = "\n".join(
            f"{i+1}. {pl.solu_name}:\n{pl.content}\n\n" for i,pl in enumerate(state["solution_plans"])
            )
        add_first_sys_msg = (
            "[Task] For a specific problem description and the associated candidate solving algorithms, select an appropriate algorithm and implement it in Python to resolve the problem.\n"       
            "[Goal] Produce a complete, well-structured, and executable Python script that correctly realizes the selected algorithm, with clear modular organization and quantitative outputs suitable for verification."     
            )
        human_msg_coder += (
                f"[Expanded Problem]:\n{state['expanded_prob']}\n\n"
                f"[Solution Plans]:\n{plans_text}"
                )
    else:
        add_first_sys_msg = (
            "[Task] For a specific problem description, implement an appropriate algorithm in Python to resolve the problem.\n"       
            "[Goal] Produce a complete, well-structured, and executable Python script that correctly realizes the solution, with clear modular organization and quantitative outputs suitable for verification."     
            )


    add_rev_sys_msg = (
                    "[Task] Revise the [Previous Code] based on the provided execution results and structured review comments.\n"       
                    "[Goal] When revising the [Previous Code], you MUST:\n"
                    "- Analyze the error trace to identify the underlying root cause of the previous failure.\n"
                    "- Address Recommendations from the [Review Comments] EXACTLY."
                    )


    ###     if Revised, after fist time
    is_revision = state.get("iteration_count", 0) > 0
    if is_revision:
        last_code = state["python_codes"][-1] if state["python_codes"] else "N/A"
        last_Runtime_opt = state["runtime_outputs"][-1] if state["runtime_outputs"] else "N/A"
        last_comments = state["review_comments"][-1] if state["review_comments"] else "N/A"
        last_tech_details = state["technical_spec"][-1] if state["technical_spec"] else "N/A"

        feedback_text = (
            f"[Technical Details]:\n{last_tech_details}\n"
            f"[Previous Code]:\n{last_code}\n"
            f"[Runtime Output]:\n{last_Runtime_opt}\n"
            f"[Review Comments] need addressed:\n{last_comments}"
        )
        config_dict = {
            "llm_name": PROG_name,
            "system_content": sys_msg_coder + add_rev_sys_msg,
            "human_content": (
                f"[Problem Statement]:\n{state['prob_todo']}\n"            
                f"{feedback_text}"
            ),
        }
        progm_response = llm_call_(config_dict, Progarmmer_module)


        # display(Markdown(f"#### 📝 Python script revised:\n{progm_response.python_codes}"))
        collapsible_mardown(f"Python script revised{state["iteration_count"]} 📝:", progm_response.python_codes)
        return {
            "technical_spec": [progm_response.technical_spec],
            "python_codes": [progm_response.python_codes],}


    ###     First time
    config_dict = {
        "llm_name": PROG_name,
        "system_content": sys_msg_coder + add_first_sys_msg,
        "human_content": human_msg_coder,
    }
    progm_response = llm_call_(config_dict, Progarmmer_module)

    # display(Markdown(f"#### 🗒️ Python script:\n{progm_response.python_codes}"))
    collapsible_mardown("Python script 🗒️:", progm_response.python_codes)
    return {
        "technical_spec": [progm_response.technical_spec],
        "python_codes": [progm_response.python_codes],}


def run_codes(state: State):
    display(Markdown("#### ⚙️ Executing"))
    Last_code_to_run = state["python_codes"][-1]  # 运行最新的代码
    runtime_output = execute_code_tool(Last_code_to_run)
    # display(Markdown(f"#### 🖥️ Runtime outputs:"))
    # print(f"{runtime_output}")
    collapsible_mardown("Runtime outputs 🖥️:", runtime_output)

    #Related images: show and review
    latest_img = load_image_tool()

    img_review_text = ""
    content_img_prompt = []
    if MULTI_MODAL_REVIEW and latest_img:
        image_rev_sys_msg = (
            "You are a Scientific Visualization Expert and Numerical Analyst.\n"
            "Analyze the generated images based on the problem statement and the corresponding Python script."
            )
        intro_msg = (
            f"[Problem statement]:\n{state['prob_todo']}\n"
            f"[Submitted Code]:\n{Last_code_to_run}\n"
            f"[Runtime Output]:\n{runtime_output}\n"
            f"Review the following {len(latest_img)} provided images\n"
        )
        content_img_prompt.append({"type": "text", "text": intro_msg})

        for n,img_path in enumerate(latest_img):
            encoded_img = encode_image_tool(img_path)
            file_name = os.path.basename(img_path)
            content_img_prompt.append({"type":"text", "text": f"Image {n+1} Filename: {file_name}\n"})
            content_img_prompt.append(
                {
                "type": "image_url",
                "image_url": {
                    "url": f"data:image/jpeg;base64,{encoded_img}",
                    }
                })
            
        display(Markdown("### 🔭🗂️ Multimodal Image Review:"))
        config_dict = {
            "llm_name": "gpt-4.1-mini",
            "system_content": image_rev_sys_msg,
            "human_content": content_img_prompt,
        }
        response_img_review = llm_call_(config_dict, Reivew_multimodal_image)


        for m,desc in enumerate(response_img_review.image_description):
                img_review_text += f"Image {m+1} Filename: {desc.img_name}\n"
                img_review_text += f"\n[Description]: {desc.img_content}\n"
                
                ### Issues list
                if desc.img_issue:
                    issues_list = []
                    for i, issue in enumerate(desc.img_issue):
                        
                        issues_list.append(f"{i+1}. {issue.severity}\n{issue.issue_item}")
                    img_review_text += f"\n[Identified Issues]:\n" + "\n".join(issues_list) + "\n\n"
                else:
                    img_review_text += "Identified Issues: None\n\n"

        # display(Markdown(f"#### 🖼️ Image-based Review:\n{img_review_text}"))
        collapsible_mardown("Image-based Review 🖼️:", img_review_text)

    return {
        "runtime_outputs": [runtime_output],
        "rev_image_description": [img_review_text]
    }


def llm_call_r3(state: State):
    display(Markdown("### 🧐 Reviewer:"))
    latest_code = state["python_codes"][-1]
    latest_output = state["runtime_outputs"][-1]
    latest_techSpec = state["technical_spec"][-1]
    iter_count = state.get("iteration_count", 0) + 1
    
    latest_img_review = state["rev_image_description"][-1] if state.get("rev_image_description") and MULTI_MODAL_REVIEW else ""

    sys_msg_reviewer = (
        "You are a Code Reviewer and Scientific Computing Expert.\n"
        "[Task] Review the reliability of numerical results and the quality of code implementation.\n"
        "[Goal] Conduct a structured review assessment based on the problem description, submitted code, and runtime output.\n"
        "Be pragmatic in your decision: if results are reasonable, accept it. But if you request a REVISE, you must be thorough and list all technical blockers."
    )

    content_parts = [f"[Problem Statement]:\n{state['prob_todo']}\n"]
    if ConsModule_ON:
        content_parts.append(f"[Expanded Problem]:\n{state['expanded_prob']}\n")
    content_parts.extend([
        f"[Code Submitted]:\n{latest_code}\n",
        f"[Technical Details]:\n{latest_techSpec}\n",
        f"[Runtime Output]:\n{latest_output}\n",
        f"[Image-based Review]:\n{latest_img_review}"
    ])
    human_msg_content = "".join(content_parts)

    config_dict = {
        "llm_name": Rev_name,
        "system_content": sys_msg_reviewer,
        "human_content": human_msg_content,
    }
    review_response = llm_call_(config_dict, Reviewer_module)


    comments_part_ = "\n".join(
        f"{i+1}. {comment.severity}\n\nCategory: {comment.category}\nIssue: {comment.issue}\n\nRecommendation: {comment.recommendation}\n\n"
        for i, comment in enumerate(review_response.review_comments)
    )
    comments_text = (
        f"Summary:\n{review_response.review_summary}\n\n"
        f"Comments:\n{comments_part_}"
    )

    display(Markdown(f"#### 📧📨📥 Decision: {review_response.review_decision.upper()}"))
    # display(Markdown(f"#### 📝 Review Comments:\n{comments_text}"))
    collapsible_mardown("Review Details 📝:", comments_text)

    return {
        "review_decision": review_response.review_decision,
        "review_comments": [comments_text],
        "iteration_count": iter_count}


def decide_next_step(state: State):
    decision = state["review_decision"]
    iter_count = state.get("iteration_count", 0)

    if decision == "accept":
        display(Markdown(f"#### ✅ Accepted after {iter_count}/{MAX_ITERATIONS} iters."))
        return "Accepted"

    if iter_count > MAX_ITERATIONS:
        display(Markdown(f"#### 🛑 Max iterations {iter_count}/{MAX_ITERATIONS} reached. Stopping."))
        return "Accepted"

    display(Markdown(f"#### 🔄 Revision requested (Iteration {iter_count}/{MAX_ITERATIONS}). Back to Programmer."))
    return "Revised"