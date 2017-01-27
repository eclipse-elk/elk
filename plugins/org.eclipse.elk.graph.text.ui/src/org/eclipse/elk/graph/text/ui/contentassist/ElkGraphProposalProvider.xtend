/*******************************************************************************
 * Copyright (c) 2016 TypeFox GmbH (http://www.typefox.io) and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *******************************************************************************/
package org.eclipse.elk.graph.text.ui.contentassist

import com.google.common.base.Predicate
import com.google.inject.Inject
import com.google.inject.Provider
import java.util.List
import org.eclipse.elk.core.data.ILayoutMetaData
import org.eclipse.elk.core.data.LayoutAlgorithmData
import org.eclipse.elk.core.data.LayoutMetaDataService
import org.eclipse.elk.core.data.LayoutOptionData
import org.eclipse.elk.core.options.CoreOptions
import org.eclipse.elk.graph.ElkEdge
import org.eclipse.elk.graph.ElkEdgeSection
import org.eclipse.elk.graph.ElkGraphElement
import org.eclipse.elk.graph.ElkLabel
import org.eclipse.elk.graph.ElkNode
import org.eclipse.elk.graph.ElkPort
import org.eclipse.elk.graph.impl.ElkPropertyToValueMapEntryImpl
import org.eclipse.elk.graph.text.services.ElkGraphGrammarAccess
import org.eclipse.emf.ecore.EObject
import org.eclipse.jface.viewers.StyledString
import org.eclipse.xtend.lib.annotations.FinalFieldsConstructor
import org.eclipse.xtext.Assignment
import org.eclipse.xtext.CrossReference
import org.eclipse.xtext.Keyword
import org.eclipse.xtext.RuleCall
import org.eclipse.xtext.conversion.impl.IDValueConverter
import org.eclipse.xtext.resource.IEObjectDescription
import org.eclipse.xtext.ui.IImageHelper
import org.eclipse.xtext.ui.editor.contentassist.ContentAssistContext
import org.eclipse.xtext.ui.editor.contentassist.ICompletionProposalAcceptor
import org.eclipse.xtext.util.Strings

import static extension org.eclipse.emf.ecore.util.EcoreUtil.*
import static extension org.eclipse.xtext.EcoreUtil2.*
import org.eclipse.elk.graph.util.ElkGraphUtil
import org.eclipse.elk.core.data.LayoutOptionData.Type

/**
 * Special content assist proposals for the ELK Graph language.
 */
class ElkGraphProposalProvider extends AbstractElkGraphProposalProvider {
    
    static val DISABLED_KEYWORDS = #{'}', ']'}
    
    @Inject IImageHelper imageHelper
    
    ElkGraphGrammarAccess grammar
    
    IDValueConverter idValueConverter
    
    @Inject
    def void initialize(Provider<IDValueConverter> idValueConverterProvider, ElkGraphGrammarAccess grammarAccess) {
        this.idValueConverter = idValueConverterProvider.get => [
            rule = grammarAccess.IDRule
        ]
        this.grammar = grammarAccess
    }
    
    override completeKeyword(Keyword keyword, ContentAssistContext context, ICompletionProposalAcceptor acceptor) {
        if (!DISABLED_KEYWORDS.contains(keyword.value) && keyword.value != context.prefix)
            super.completeKeyword(keyword, context, acceptor)
    }
    
    override protected doCreateStringProposals() {
        false
    }

    override complete_PropertyKey(EObject model, RuleCall ruleCall, ContentAssistContext context,
            ICompletionProposalAcceptor acceptor) {
        switch model {
            ElkNode: {
                if (model.parent === null || !model.children.empty)
                    proposeProperties(model, model.algorithm, LayoutOptionData.Target.PARENTS, context, acceptor)
                if (model.parent !== null)
                    proposeProperties(model, model.parent.algorithm, LayoutOptionData.Target.NODES, context, acceptor)
            }
            ElkEdge: {
                proposeProperties(model, model.algorithm, LayoutOptionData.Target.EDGES, context, acceptor)
            }
            ElkPort: {
                proposeProperties(model, model.algorithm, LayoutOptionData.Target.PORTS, context, acceptor)
            }
            ElkLabel: {
                proposeProperties(model, model.algorithm, LayoutOptionData.Target.LABELS, context, acceptor)
            }
        }
    }
    
    private def getAlgorithm(ElkGraphElement element) {
        var node = element.getContainerOfType(ElkNode)
        if (node !== null) {
            if ((element instanceof ElkLabel || element instanceof ElkPort) && node.parent !== null)
                node = node.parent
            val algorithmId = node.getProperty(CoreOptions.ALGORITHM)
            if (!algorithmId.nullOrEmpty)
                return LayoutMetaDataService.instance.getAlgorithmDataBySuffix(algorithmId)
        }
    }
    
    protected def proposeProperties(ElkGraphElement element, LayoutAlgorithmData algorithmData,
            LayoutOptionData.Target targetType, ContentAssistContext context, ICompletionProposalAcceptor acceptor) {
        val metaDataService = LayoutMetaDataService.instance
        val filteredOptions = metaDataService.optionData.filter[ o |
            targetType === null || o.targets.contains(targetType)
        ].filter[ o |
            algorithmData === null || algorithmData.knowsOption(o) || CoreOptions.ALGORITHM == o
        ].filter[ o |
            element === null || !element.properties.map.containsKey(o)
        ]
        for (option : filteredOptions) {
            val idSplit = Strings.split(option.id, '.')
            val prefixSplit = Strings.split(context.prefix, '.')
            var foundMatch = false
            var i = idSplit.size - 1
            if (i >= 1 && option.group == idSplit.get(i - 1)) {
                i--
            }
            while (i >= 0 && !foundMatch) {
                val suffix = idSplit.drop(i)
                if (metaDataService.getOptionDataBySuffix(suffix.join('.')) !== null && suffix.startsWith(prefixSplit))
                    foundMatch = true
                else
                    i--
            }
            if (foundMatch) {
                val suffix = idSplit.drop(i)
                val proposal = createCompletionProposal(suffix.convert, option.getDisplayString(suffix), getImage(option, null), context)
                acceptor.accept(proposal)
            }
        }
    }
    
    override completeProperty_Value(EObject model, Assignment assignment, ContentAssistContext context,
            ICompletionProposalAcceptor acceptor) {
        if (model instanceof ElkPropertyToValueMapEntryImpl) {
            val property = model.key
            if (property instanceof LayoutOptionData) {
                if (CoreOptions.ALGORITHM == property)
                    proposeAlgorithms(context, acceptor)
                else 
                    typeAwarePropertyValueProposal(property, assignment, context, acceptor)
            }
        }
    }
    
    private def typeAwarePropertyValueProposal(LayoutOptionData property, Assignment assignment, ContentAssistContext context, 
         ICompletionProposalAcceptor acceptor) {
         
         switch (property.type) {
             case Type.BOOLEAN,
             case Type.ENUM, 
             case Type.ENUMSET: {
                 val choices = property.choices
                 for (var i = 0; i < choices.length; i++) {
                    val proposal = choices.get(i)
                    val enumVal = property.getEnumValue(i)
                    
                    val displayString = new StyledString(proposal)
                    var priority = 3
                    if (ElkGraphUtil.isExperimentalPropertyValue(enumVal)) {
                        displayString.append(" - Experimental", StyledString.COUNTER_STYLER);
                        priority = 1
                    } else if (ElkGraphUtil.isAdvancedPropertyValue(enumVal)) {
                        displayString.append(" - Advanced", StyledString.COUNTER_STYLER);
                        priority = 2
                    }
                    acceptor.accept(createCompletionProposal(proposal, displayString, getImage(property, proposal), priority, "", context))
                 }
             }
             case DOUBLE:
                acceptor.accept(createCompletionProposal("0.0", property.getType().toString(), null, context))
             case INT:
                acceptor.accept(createCompletionProposal("0", property.getType().toString(), null, context))
             case OBJECT: {
                val proposal = try {
                    "\"" + property.getOptionClass().newInstance().toString() + "\"";
                } catch (InstantiationException e) ""
                  catch (IllegalAccessException e) ""
                acceptor.accept(createCompletionProposal(proposal, property.getType().toString(), null, context)) 
             }
             default: { } // nothing to propose
         }
    }
    
    protected def proposeAlgorithms(ContentAssistContext context, ICompletionProposalAcceptor acceptor) {
        val metaDataService = LayoutMetaDataService.instance
        for (algorithm : metaDataService.algorithmData) {
            val idSplit = Strings.split(algorithm.id, '.')
            val prefixSplit = Strings.split(context.prefix, '.')
            var foundMatch = false
            var i = idSplit.size - 1
            while (i >= 0 && !foundMatch) {
                val suffix = idSplit.drop(i)
                if (metaDataService.getAlgorithmDataBySuffix(suffix.join('.')) !== null && suffix.startsWith(prefixSplit))
                    foundMatch = true
                else
                    i--
            }
            if (foundMatch) {
                val suffix = idSplit.drop(i)
                val proposal = createCompletionProposal(suffix.convert, algorithm.getDisplayString(suffix), null, context)
                acceptor.accept(proposal)
            }
        }
    }
    
    private def startsWith(Iterable<String> strings, List<String> prefix) {
        if (prefix.empty)
            return true
        val stringList = strings.toList
        for (var i = 0; i < stringList.size - prefix.size + 1; i++) {
            var j = 0
            var matches = true
            while (j < prefix.size && matches) {
                matches = stringList.get(i + j).startsWith(prefix.get(j))
                j++
            }
            if (matches)
                return true
        }
        return false
    }
    
    private def convert(Iterable<String> suffix) {
        suffix.map[idValueConverter.toString(it)].join('.')
    }
    
    private def getDisplayString(ILayoutMetaData data, Iterable<String> suffix) {
        new StyledString(suffix.join('.'))
            + new StyledString(''' «'\u2013'» «data.name» («data.id»)''', StyledString.QUALIFIER_STYLER)
    }
    
    override protected getImage(EObject eObject) {
        if (eObject instanceof Keyword) {
            val key = switch eObject {
                case grammar.rootNodeAccess.graphKeyword_1_0: 'elkgraph'
                case grammar.elkNodeAccess.nodeKeyword_0: 'elknode'
                case grammar.elkEdgeAccess.edgeKeyword_0: 'elkedge'
                case grammar.elkPortAccess.portKeyword_0: 'elkport'
                case grammar.elkLabelAccess.labelKeyword_0: 'elklabel'
            }
            if (key !== null)
                return imageHelper.getImage(key + '.gif')
        }
        return super.getImage(eObject)
    }
    
    private def getImage(LayoutOptionData option, String value) {
        val key = switch option.type {
            case BOOLEAN:
                if (value == 'false') 'prop_false'
                else 'prop_true'
            case INT: 'prop_int'
            case DOUBLE: 'prop_double'
            case ENUM, case ENUMSET: 'prop_choice'
            default: 'prop_text'
        }
        return imageHelper.getImage(key + '.gif')
    }
    
    private def +(StyledString s1, StyledString s2) {
        s1.append(s2)
    }
    
    override completeElkEdgeSection_IncomingShape(EObject model, Assignment assignment, ContentAssistContext context,
            ICompletionProposalAcceptor acceptor) {
        if (model instanceof ElkEdgeSection)
            lookupCrossReference((assignment.terminal as CrossReference), context, acceptor,
                    new SectionShapeFilter(model, SectionShapeFilter.INCOMING))
        else
            super.completeElkEdgeSection_IncomingShape(model, assignment, context, acceptor)
    }
    
    override completeElkEdgeSection_OutgoingShape(EObject model, Assignment assignment, ContentAssistContext context, ICompletionProposalAcceptor acceptor) {
        if (model instanceof ElkEdgeSection)
            lookupCrossReference((assignment.terminal as CrossReference), context, acceptor,
                    new SectionShapeFilter(model, SectionShapeFilter.OUTGOING))
        else
            super.completeElkEdgeSection_OutgoingShape(model, assignment, context, acceptor)
    }
    
    override completeElkSingleEdgeSection_IncomingShape(EObject model, Assignment assignment, ContentAssistContext context, ICompletionProposalAcceptor acceptor) {
        if (model instanceof ElkEdgeSection)
            lookupCrossReference((assignment.terminal as CrossReference), context, acceptor,
                    new SectionShapeFilter(model, SectionShapeFilter.INCOMING))
        else
            super.completeElkSingleEdgeSection_IncomingShape(model, assignment, context, acceptor)
    }
    
    override completeElkSingleEdgeSection_OutgoingShape(EObject model, Assignment assignment, ContentAssistContext context, ICompletionProposalAcceptor acceptor) {
        if (model instanceof ElkEdgeSection)
            lookupCrossReference((assignment.terminal as CrossReference), context, acceptor,
                    new SectionShapeFilter(model, SectionShapeFilter.OUTGOING))
        else
            super.completeElkSingleEdgeSection_OutgoingShape(model, assignment, context, acceptor)
    }
    
    @FinalFieldsConstructor
    static class SectionShapeFilter implements Predicate<IEObjectDescription> {
        
        static val INCOMING = 0
        static val OUTGOING = 1
        
        val ElkEdgeSection section
        val int type
        
        override apply(IEObjectDescription input) {
            switch type {
                case INCOMING:
                    input.isInList(section.parent.sources)
                case OUTGOING:
                    input.isInList(section.parent.targets)
                default: true
            }
        }
        
        private def isInList(IEObjectDescription input, List<? extends EObject> list) {
            val object = input.EObjectOrProxy
            if (object.eIsProxy) {
                list.exists[URI == input.EObjectURI]
            } else {
                list.contains(object)
            }
        }
        
    }
    
}
