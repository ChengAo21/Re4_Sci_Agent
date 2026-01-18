from pydantic import BaseModel, Field
from typing_extensions import TypedDict, Literal
from typing import Annotated, List
import operator
from config import ALTERNATIVE_NUM, MAX_COMMENTS

class State(TypedDict):
    # --- Problem ---
    prob_todo: str

    # --- Design Section ---
    expanded_prob: str
    solution_plans: List[str]
    
    # --- Execution Section ---
    technical_spec: Annotated[List[str], operator.add]
    python_codes: Annotated[List[str], operator.add]
    runtime_outputs: Annotated[List[str], operator.add]
    
    # --- Review Section ---
    review_decision: Literal["accept", "revise"]
    review_comments: Annotated[List[str], operator.add]
    iteration_count: int

    # --- Review multimodal IMAGE ---
    rev_image_description: Annotated[List[str], operator.add]  


######## 
class solu_plan(BaseModel):
    solu_name: str = Field(description="Short title of the solution strategy")
    content: str = Field(description=
            "Structured formulation in concise pseudocode or numbered steps style.\n"
            "Requirements:\n"
            "- governing idea\n"
            "- algorithmic steps\n"
            "- stability / accuracy / complexity / efficiency limitations\n"
            "STRICT CONSTRAINT: \n"
            "-Do not directly provide Python script here."
            "-Stop immediately after the limitations sections."
        )

class Consultant_module(BaseModel):
    expanded_prob: str = Field( description=
        "Elaborate on the original problem outline while preserving its core meaning.\n"
        "Explicitly indentify the primary mathematical and numerical challenges."
            )
    solution_plans: List[solu_plan] = Field( description=
        f"Provide {ALTERNATIVE_NUM} alternative solution strategies. Do not write long continuous paragraphs."
        )



########
class Progarmmer_module(BaseModel):
    technical_spec: str = Field(description = 
            "Concise architectural overview for rapid comprehension of core algorithm and data flow."
    )
    python_codes: str = Field(description=
            "A complete, bug-free Python script.\n"
            "Requirements:\n"
            "- Enclosed in a single ```python ... ``` markdown block.\n"
            "- The content INSIDE the markdown block MUST be PURE PYTHON CODE ONLY.\n"
            "- All required imports included at the top.\n"            
            "- Modular structure with clear function definitions.\n"
            "- Use brief, clear code comments, avoid verbose or tutorial-style explanations.\n"
            "- Do NOT rebuild or modify matrices inside the time loop; do it ONCE at the start.\n"
            "- Must include quantitative outputs: printed metrics (e.g., error norm) or labeled plots (with axes, title)."
        )


########
class review_spec(BaseModel):
    category: Literal["runtime","structure","accuracy",] = Field( description=
            "Specify the type of issue being addressed.\n"
            "- runtime: Execution errors or failures\n"
            "- structure: Code organization or quality\n" 
            "- accuracy: Numerical stability, precision, convergence, or error behavior\n"
            "MAJOR/MINOR Issues inferred from image-based review MUST be re-categorized based on their underlying cause (runtime, structure, or accuracy)."
            )
    issue: str = Field( description= "A concrete issue identified.")
    severity: Literal["MINOR", "MAJOR"] = Field( description=
            "Severity of the issue:\n"
            "- MINOR: Cosmetic, stylistic, or secondary issues that do not affect correctness or results\n"
            "- MAJOR: Issues that cause runtime failure, non-physical behavior, incorrect conclusions, or invalid results."
            )
    recommendation: str = Field( description= "Specific suggestions for improvement and refinement.")

class Reviewer_module(BaseModel):
    review_summary: str = Field( description=
            "A concise summary determining whether the employed algorithm is appropriate and if the problem is perfectly solved."
            ) 
    review_decision: Literal["accept", "revise"] = Field( description=
            "Accept only if:\n"         
            "- code runs without errors\n"
            "- numerical results are correct and quantitatively accurate\n"
            "- no MAJOR issues are identified."
            "Revise only if:\n"
            "- there is any runtime error or the results are fundamentally wrong.\n"
            "- any MAJOR issue exists, including those inferred from image-based review."
            )
    review_comments: List[review_spec] = Field( description=
            f"Provide up to {MAX_COMMENTS} specific review comments.\n"
            "If accepting, leave list empty or provide optional suggestions.\n"
            "If revising, provide a comprehensive list covering all identified issues.\n"
            "Each comment must identify a concrete issue.\n"
            "Vague or ambiguous feedback is not allowed."
            )



######## 
class ImageIssueItem(BaseModel):
    issue_item: str = Field(description=
        "Concrete problem or discrepancy observed in the image.\n"
        "MUST refer to incorrect or implausible trends, magnitudes, spatial patterns, continuity, physical behaviors."
    )
    severity: Literal["MINOR", "MAJOR"] = Field(description=
        "Severity of the issue:\n"
        "- MINOR: Cosmetic or secondary issue. Does not affect physical validity or main results.\n"
        "- MAJOR: Non-physical behavior or clearly abnormal image indicating numerical or modeling failure."
    )

class single_rev_img_spec(BaseModel):
    img_name: str = Field(description=
        "Filename of the reviewed image."
    )
    img_content: str = Field(description=
        "Evaluation of the image, focus on:\n"
        "- consistency with problem requirements\n"
        "- primary quantities or fields shown"
    )
    img_issue: List[ImageIssueItem] = Field(description=
        "List of issues. Leave empty if no concrete problems are found."
    )

class Reivew_multimodal_image(BaseModel):
    image_description: List[single_rev_img_spec] = Field(description=
        "Structured description of each provided image."
    )