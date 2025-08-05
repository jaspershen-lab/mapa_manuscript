## Task Description
You will be given a meta-module containing multiple functional modules from different tissues. Your goal is to:

1. **Generate a concise, biologically meaningful meta-module name** (3-8 words)
2. **Provide a comprehensive summary** explaining the biological significance of this meta-module.

## Input Information Format
For each functional module in the meta-module, you will receive:
- **Module Summary**: Functional description of the module
- **Tissue**: The tissue where this module was identified

## Analysis Guidelines

### For Meta-module Naming:
- Focus on the **core biological process** or **pathway theme** shared across functional modules
- Consider **tissue specificity** if relevant (e.g., "Cardiac-specific", "Multi-tissue")
- Use standard biological terminology
- Keep it concise but informative

### For Summary Generation:
- **Identify common biological themes** across all modules in the meta-module
- **Mention key pathways or processes** represented in the meta-module
- **Consider cross-tissue communication** or shared regulatory mechanisms
- **Highlight their biological relevance**

## Output Format
Respond with a valid JSON object using exactly this format (no markdown, no code blocks, no additional text):

{"meta-module_name": "Your concise meta-module name here", "summary": "A comprehensive biological summary explaining the meta-module's significance, tissue patterns, and potential biological implications"}

Start your response directly with the opening brace { and end with the closing brace }.

## **Actual Input for Generation**
**Functional modules information in this meta-modules**:
{module_data}

## Example Analysis Approach
1. Read through all functional module summaries to identify common themes
2. Note tissue distribution
3. Consider biological relationships between tissues
4. Formulate a name that captures the essence of the meta-module
5. Write a summary that explains the biological story this meta-module tells

Please analyze the provided meta-modules and generate the meta-module name and summary following the guidelines above.
