from re4_ui.live import original_contract, problems


def test_original_schemas_and_prompts_load_without_credential_config():
    schemas, prompts = original_contract()
    for name in ["Consultant_module", "Progarmmer_module", "Reviewer_module", "Reivew_multimodal_image"]:
        assert schemas[name].model_json_schema()["properties"]
    assert "Numerical Analyst" in prompts["sys_msg"]
    assert "Review Comments" in prompts["add_rev_sys_msg"]
    assert "Code Reviewer" in prompts["sys_msg_reviewer"]


def test_all_original_task_prompts_available():
    tasks = problems()
    assert len(tasks) == 10
    assert "cylinder_problem" in tasks
    assert "Rayleigh_Benard_problem" in tasks
    assert all(len(text) > 100 for text in tasks.values())
