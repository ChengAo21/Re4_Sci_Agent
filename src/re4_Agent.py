import utils
from prompts_template import prompts_roles
from llm_agent import LLMAgent
from code_executor import CodeProcessor
from FileHandler import File_Handler



class re4_SciComp_Agent:
    def __init__(self,
                 original_problem: str,
                 consultant_config: dict = None,
                 programmer_config: dict = None,
                 reviewer_config: dict = None,
                 time_rev: int = 2,
                 path_save: str = "./results_re4/",                 
                 ):
        
        self.original_problem = original_problem
        self.time_review = time_rev
        
        self.consutant_agent = self._init_llm_agent(consultant_config,)
        self.programmer_agent = self._init_llm_agent(programmer_config,)
        self.reviewer_agent = self._init_llm_agent(reviewer_config,)

        self.code_executor = CodeProcessor()
        self.file_handler = File_Handler(path_save)

    def _init_llm_agent(self, model_config: dict,):
        agent = LLMAgent(model_config['model'], model_config['api_key'], model_config.get('base_url', None))
        print(f"[InIt Roles] {model_config['role']} with LLM: {model_config['model']}\n")
        return agent
    
    def step1_augument_problem(self,have_done = False) -> str:
        print(f"[Step 1]: Consultant expand the problem. Just read? {have_done}\n")
        if not have_done:
            prompt_forConsultant = prompts_roles.generate_prompt_for_Consultant(self.original_problem,)
            expanded_problem  = self.consutant_agent.execute(prompt_forConsultant,)
            self.file_handler.write_content(expanded_problem, "1_Consultant_augumented_problem", 0)
        else:
            expanded_problem = self.file_handler.read_content("1_Consultant_augumented_problem", 0)

        print(f"[re4_Agent] Consultant's agumentation done. Just read? {have_done}\n")
        return expanded_problem
    
    def step2_initial_solution(self, problem_analysis: str) -> str:
        print(f'[Step 2]: Programmer generate initial solution\n')
        prompt_forInit_Programmer = prompts_roles.generate_prompt_for_Programmer_Initial(self.original_problem, problem_analysis)
        initial_solution = self.programmer_agent.execute(prompt_forInit_Programmer,)

        exe_result = self.code_executor.run_codes(initial_solution)
        self.file_handler.write_content_code_out(exe_result, "2_Programmer_initial_solution", 0)
        print(f"[re4_Agent] Programmer's initial solution done.\n")
        return exe_result
    
    def step3_review_answer(self, problem_analysis: str, results_execution: str, count_times: int = 0) -> str:
        print(f"[Step 3]: Reviewer review the runtime outputs\n")
        prompt_forReviewer = prompts_roles.generate_prompt_for_Reviewer(self.original_problem, problem_analysis, results_execution)
        review_feedback = self.reviewer_agent.execute(prompt_forReviewer,)
        self.file_handler.write_content(review_feedback, "3_Reviewer_feedback", count_times)
        print(f"[re4_Agent] Reviewer's feedback done.\n")
        return review_feedback
    
    def step4_revise_solution(self, previous_solu: str, review_feedback: str, count_times: int = 0) -> str:
        print(f"[Step 4]: Programmer revise the solution based on reviewer's feedback\n")
        prompt_forRevise_Programmer = prompts_roles.generate_prompt_for_Programmer_Revision(self.original_problem, previous_solu, review_feedback)
        revised_solution = self.programmer_agent.execute(prompt_forRevise_Programmer,)

        exe_result = self.code_executor.run_codes(revised_solution)
        self.file_handler.write_content_code_out(exe_result, "4_Programmer_revised_solution", count_times)
        print(f"[re4_Agent] Programmer's revised solution done.\n")
        return exe_result
    
    def run(self) -> None:

        print("="*50)
        print("🤖: re4 Agent Launch🚀")
        print("="*50+'\n')

        # 1. 扩展问题
        agumented_problem = self.step1_augument_problem(have_done=True)

        # 2. 生成初始代码
        programmer_results = self.step2_initial_solution(agumented_problem)

        # 3. 迭代评审与修订
        for i in range(1,self.time_review+1):
            review_feedback = self.step3_review_answer(agumented_problem, programmer_results, i)
            revised_results = self.step4_revise_solution(programmer_results, review_feedback, i)
            programmer_results = revised_results


        print("\n" + "="*50)
        print("✅ All steps completed. Results saved.")
        print(f"✍️Next steps are to release the Streamlit version and apply LangChain.😎🤓")
        print("="*50)    

