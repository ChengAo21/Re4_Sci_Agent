# Re4-Scientific-Computing-Agent-with-Rewriting-Resolution-Review-and-Revision
Large language models (LLMs) serve as an active and promising field of generative artificial intelligence and have demonstrated abilities to perform complex tasks in multiple domains, including mathematical and scientific reasoning. In this work, we construct a novel agent framework for solving representative problems in scientific computing. The proposed agent, incorporating a "rewriting-resolution-review-revision" logical chain via three reasoning LLMs (functioning as the Consultant, Reviewer, and Programmer, respectively), is integrated in a collaborative and interactive manner. The Consultant module endows the agent with knowledge transfer capabilities to link problems to professional domain insights, thereby rewriting problem descriptions through text augmentation. The Programmer module is responsible for generating and executing well-structured code to deliver the problem resolution. The Reviewer module equips the agent with the capacity for self-debugging and self-refinement through interactive feedback with code runtime outputs. By leveraging the end-to-end review mechanism, the executable code provided by the Programmer attains the iterative revision. A comprehensive evaluation is conducted on the performance of the proposed agent framework in solving PDEs, ill-conditioned linear systems, and data-driven physical analysis problems. Compared to single-model, this collaborative framework significantly improves the bug-free code generation rate and reduces the occurrence of non-physical solutions, thereby establishing a highly reliable framework for autonomous code generation based on natural language descriptions. The review mechanism improved the average execution success (bug-free code and non-NaN solutions) rate of the latest reasoning models. In summary, our agent framework establishes automatic code generation and review as a promising scientific computing paradigm.

If you find this repository useful, please consider giving a star ⭐ and cite our paper.
```
@article{cheng2025re4,
  title={Re4: Scientific Computing Agent with Rewriting, Resolution, Review and Revision},
  author={Cheng, Ao and Zhang, Lei and He, Guowei},
  journal={arXiv preprint arXiv:2508.20729},
  year={2025}
}
```
