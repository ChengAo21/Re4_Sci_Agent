class prompts_roles:

    @staticmethod
    def generate_prompt_for_Consultant(original_problem: str) -> str:
        
        task_description = "Expand the context of a scientific problem and provide a variety of solution options.".strip()
        
        goal_description = "For specific problems, expand the problem background while maintaining the original meaning and carefully analyse the problems, then provide multiple detailed solution algorithm.".strip()

        output_requirement = """
        Provide a robust context of the problem, maintaining the same main context, and offer multiple detailed solution plans.
        1. Expand the original problem's context as much as possible based on its simple outline.
        2. Please provide a variety of detailed solution plans, with the requirement that each plan is thorough and specific.
        """.strip()

        prompt_ = f"""
        [Task]: {task_description}
        [Goal]: {goal_description}

        Original Specific Problem: {original_problem}

        [Output Requirement]: {output_requirement}
        """.strip()

        return prompt_

    @staticmethod
    def generate_prompt_for_Programmer_Initial(original_problem: str, problem_analysis: str) -> str:

        task_description = "For the specific problem, select an appropriate solving algorithm and generate Python code to complete the problem-solving process.".strip()
        goal_description = "Generate structured Python code by first leveraging your understanding based on the detailed descriptions of the specific problem and its solving algorithm.".strip()

        output_requirement = """
        1. Based on the specific description of the problem and solving algorithms, select an appropriate framework to accomplish the problem-solving.
        2. Ensure the code includes:
        - Ensure only use # for annotations to explain each section of the code.
        - Technical explanation for the solving algorithm.
        - Implement detailed result printing within the code architecture.
        3. The output code should be enclosed with ``` as follows:
        ```python
        # Technical explanation content

        #Hints:
        #Properly handle the dimensions, indexing, and related operations of arrays or tensors.
        #Provide a bug-free and well-structured code for reproductibility.
        **Python code here**

        # Optimized explanation content
        ```
        """.strip()

        prompt_ = f"""
        [Task]: {task_description}
        [Goal]: {goal_description}

        [Original Specific Problem]: {original_problem}

        [Consultant's expansion]: {problem_analysis}

        [Output Requirement] {output_requirement}
        """.strip()

        return prompt_

    @staticmethod
    def generate_prompt_for_Reviewer(original_problem: str, problem_analysis: str, results_Fromcode: str) -> str:

        task_description = "Review the answer provided by the programmer to the specific problem and provide feedback".strip()
        goal_description = "Please conduct a detailed review according to the complete description of a scientific problem and analyze the results output by the programmer's code.".strip()

        output_requirement = """
        1. Based on the detailed analysis of the problem and the Programmer’s solution, analyze the results output by the code and provide detailed feedback, guide the programmer in further deepening his understanding of the problem and solving it with greater perfection. 
        2. Ensure your feedback includes:
        - Determine whether the programmer has perfectly solved the problem and whether the most appropriate algorithm has been used.
        - Assist the programmer in checking and refining runtime errors and warnings in the code.
        - Provide suggestions for optimizing the code, including but not limited to algorithm optimization, code structure optimization, and handling of potential errors in the code.
        - Your feedback includes posteriori issue identification based on programmer's results, and may also include a priori analysis based on your understanding of the specific problem.
        """.strip()

        prompt_ = f"""
        [Task]: {task_description}
        [Goal]: {goal_description}

        [Original Specific Problem]: {original_problem}

        [Consultant's expansion]: {problem_analysis}

        [Programmer's solution]: {results_Fromcode}

        [Output Requirement] {output_requirement}
        """.strip()

        return prompt_

    @staticmethod
    def generate_prompt_for_Programmer_Revision(original_problem: str, previous_solu: str, review_feedback: str) -> str:

        task_description = "For the specific problem, select an appropriate solving algorithm and generate Python code to complete the problem-solving process.".strip()
        goal_description = "Generate structured Python code by first leveraging your understanding based on the detailed descriptions of the specific problem and its solving algorithm.".strip()

        output_requirement = """
        1. Generate refined Python code by rationally incorporating feedback from the reviewer model, based on your comprehensive understanding of the specific problem and its solving algorithms.
        2. Ensure the code includes:
        - Ensure only use # for annotations to explain each section of the code.
        - Technical explanation for the solving algorithm.
        - Implement detailed result printing within the code architecture.
        - Explain the specific parts of the code that have been optimized.
        3. The output code should be enclosed with ``` as follows:
        ```python
        # Technical explanation content

        #Hints:
        #Properly handle the dimensions, indexing, and related operations of arrays or tensors.
        #Provide a bug-free and well-structured code for reproductibility.
        **Python code here**

        # Optimized explanation content
        ```
        """.strip()

        prompt_ = f"""
        [Task]: {task_description}
        [Goal]: {goal_description} Then fully incorporating feeedback from the reviewer model

        [Original Specific Problem]: {original_problem}

        [Programmer's previous solution]: {previous_solu}

        [Reviewer's Feedback]: {review_feedback}

        [Output Requirement] {output_requirement}
        """.strip()

        return prompt_
