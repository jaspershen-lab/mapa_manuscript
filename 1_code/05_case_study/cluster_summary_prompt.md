## Task Description
You will be given a cluster containing multiple functional modules identified from bulk RNA-seq data across different tissues. Each functional module represents a group of co-regulated pathways that are either upregulated or downregulated in specific tissues. Your goal is to:

1. **Generate a concise, biologically meaningful cluster name** (3-8 words)
2. **Provide a comprehensive summary** explaining the biological significance of this cluster

## Input Information Format
For each functional module in the cluster, you will receive:
- **Module Summary**: Functional description of the module
- **Tissue**: The tissue where this module was identified
- **Regulation Direction**: "Cluster U" for "upregulated" or "Cluster D" for "downregulated" with aging

## Analysis Guidelines

### For Cluster Naming:
- Focus on the **core biological process** or **pathway theme** shared across modules
- Consider **tissue specificity** if relevant (e.g., "Cardiac-specific", "Multi-tissue")
- Include **regulation pattern** if consistent (e.g., "Upregulated", "Stress-induced")
- Use standard biological terminology
- Keep it concise but informative

### For Summary Generation:
- **Identify common biological themes** across all modules in the cluster
- **Highlight tissue-specific patterns** and their biological relevance
- **Discuss regulation patterns** (consistent up/down-regulation vs. tissue-specific differences)
- **Explain potential biological significance** (e.g., coordinated response, developmental program, disease association)
- **Consider cross-tissue communication** or shared regulatory mechanisms
- **Mention key pathways or processes** represented in the cluster

## Output Format
Respond with a valid JSON object using exactly this format (no markdown, no code blocks, no additional text):

{"module_name": "Your concise module name here", "summary": "A comprehensive biological summary explaining the cluster's significance, tissue patterns, regulation directions, and potential biological implications"}

Start your response directly with the opening brace { and end with the closing brace }.

## **Actual Input for Generation**
**Functional Modules information in this Cluster**:
{module_data}

## Example Analysis Approach
1. Read through all module summaries to identify common themes
2. Note tissue distribution and regulation patterns
3. Consider biological relationships between tissues
4. Formulate a name that captures the essence of the cluster
5. Write a summary that explains the biological story this cluster tells

Please analyze the provided cluster and generate the cluster name and summary following the guidelines above.
