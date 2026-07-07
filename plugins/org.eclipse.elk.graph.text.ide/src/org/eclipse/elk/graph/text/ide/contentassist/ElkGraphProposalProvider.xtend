/*******************************************************************************
 * Copyright (c) 2017 TypeFox GmbH (http://www.typefox.io) and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.text.ide.contentassist

import com.google.common.base.Predicate
import com.google.inject.Inject
import com.google.inject.Provider
import java.util.List
import org.eclipse.elk.core.data.ILayoutMetaData
import org.eclipse.elk.core.data.LayoutDataContentAssist
import org.eclipse.elk.core.data.LayoutOptionData
import org.eclipse.elk.core.options.CoreOptions
import org.eclipse.elk.graph.ElkEdgeSection
import org.eclipse.elk.graph.ElkGraphElement
import org.eclipse.elk.graph.impl.ElkPropertyToValueMapEntryImpl
import org.eclipse.elk.graph.text.services.ElkGraphGrammarAccess
import org.eclipse.emf.ecore.EObject
import org.eclipse.xtend.lib.annotations.FinalFieldsConstructor
import org.eclipse.xtext.Assignment
import org.eclipse.xtext.CrossReference
import org.eclipse.xtext.Keyword
import org.eclipse.xtext.conversion.impl.IDValueConverter
import org.eclipse.xtext.ide.editor.contentassist.ContentAssistContext
import org.eclipse.xtext.ide.editor.contentassist.ContentAssistEntry
import org.eclipse.xtext.ide.editor.contentassist.IIdeContentProposalAcceptor
import org.eclipse.xtext.ide.editor.contentassist.IProposalConflictHelper
import org.eclipse.xtext.ide.editor.contentassist.IdeContentProposalProvider
import org.eclipse.xtext.resource.IEObjectDescription
import org.eclipse.xtext.util.Strings

import static extension org.eclipse.emf.ecore.util.EcoreUtil.*

/**
 * Special content assist proposals for the ELK Graph language.
 */
class ElkGraphProposalProvider extends IdeContentProposalProvider {
    
    static val DISABLED_KEYWORDS = #{'}', ']'}
    
    ElkGraphGrammarAccess grammar
    
    IDValueConverter idValueConverter
    
    @Inject IProposalConflictHelper conflictHelper
    
    @Inject
    def void initialize(Provider<IDValueConverter> idValueConverterProvider, ElkGraphGrammarAccess grammarAccess) {
        this.idValueConverter = idValueConverterProvider.get => [
            rule = grammarAccess.IDRule
        ]
        this.grammar = grammarAccess
    }
    
    override protected filterKeyword(Keyword keyword, ContentAssistContext context) {
        !DISABLED_KEYWORDS.contains(keyword.value) && keyword.value != context.prefix
    }
    
    override protected _createProposals(Keyword keyword, ContentAssistContext context, IIdeContentProposalAcceptor acceptor) {
        if (filterKeyword(keyword, context)) {
            val entry = proposalCreator.createProposal(keyword.value, context)
            if (entry !== null) {
                entry.kind = ContentAssistEntry.KIND_KEYWORD
                entry.source = keyword
                acceptor.accept(entry, proposalPriorities.getKeywordPriority(keyword.value, entry))
            }
        }
    }
    
    override protected _createProposals(Assignment assignment, ContentAssistContext context, IIdeContentProposalAcceptor acceptor) {
        switch assignment {
            case grammar.propertyAccess.keyAssignment_0:
                completePropertyKey(context, acceptor)
            case grammar.propertyAccess.valueAssignment_2_0,
            case grammar.propertyAccess.valueAssignment_2_1,
            case grammar.propertyAccess.valueAssignment_2_2,
            case grammar.propertyAccess.valueAssignment_2_3:
                completePropertyValue(context, acceptor)
            default:
                super._createProposals(assignment, context, acceptor)
        }
    }
    
    protected def void completePropertyKey(ContentAssistContext context, IIdeContentProposalAcceptor acceptor) {
        if (conflictHelper.existsConflict('a', context)) {
            // Early-exit in case any property id would conflict with the previous token
            return
        }
        val model = context.currentModel
        if (model instanceof ElkGraphElement) {
            LayoutDataContentAssist.getLayoutOptionProposals(model, context.prefix).forEach [ p |
                val entry = new ContentAssistEntry => [
                    proposal = Strings.split(p.proposal, '.').convert
                    prefix = context.prefix
                    kind = ContentAssistEntry.KIND_PROPERTY
                    label = p.label ?: p.proposal
                    description = getDescription(p.data)
                    documentation = p.data.description
                    source = p.data
                ]
                acceptor.accept(entry, proposalPriorities.getDefaultPriority(entry))
            ]
        }
    }
    
    /** There's another option that allows to select a layout algorithm other than {@link CoreOptions#ALGORITHM}. 
     *  To avoid a dependency to that plugin, the option is hard-coded here. */
    public static val DISCO_LAYOUT_ALG_ID = "org.eclipse.elk.disco.componentCompaction.componentLayoutAlgorithm";
    
    protected def void completePropertyValue(ContentAssistContext context, IIdeContentProposalAcceptor acceptor) {
        val model = context.currentModel
        if (model instanceof ElkPropertyToValueMapEntryImpl) {
            val option = model.key
            if (option instanceof LayoutOptionData) {
                if (CoreOptions.ALGORITHM == option || option.id == DISCO_LAYOUT_ALG_ID)
                    proposeAlgorithms(context, acceptor)
                else
                    typeAwarePropertyValueProposal(option, context, acceptor)
            }
        }
    }
    
    private def typeAwarePropertyValueProposal(LayoutOptionData property, ContentAssistContext context,
        IIdeContentProposalAcceptor acceptor) {
        LayoutDataContentAssist.getLayoutOptionValueProposal(property, context.prefix).forEach [ p |
            val entry = proposalCreator.createProposal(p.proposal, context, ContentAssistEntry.KIND_VALUE) [
                label = p.label
                source = property
            ]
            acceptor.accept(entry, proposalPriorities.getDefaultPriority(entry))
        ]
    }
    
    protected def proposeAlgorithms(ContentAssistContext context, IIdeContentProposalAcceptor acceptor) {
        if (conflictHelper.existsConflict('a', context)) {
            // Early-exit in case any algorithm id would conflict with the previous token
            return
        }
        LayoutDataContentAssist.getLayoutAlgorithmProposals(context.prefix).forEach [ p |
            val entry = new ContentAssistEntry => [
                proposal = Strings.split(p.proposal, '.').convert
                prefix = context.prefix
                kind = ContentAssistEntry.KIND_VALUE
                label = p.label ?: p.proposal
                description = getDescription(p.data)
                documentation = p.data.description
                source = p.data
            ]
            acceptor.accept(entry, proposalPriorities.getDefaultPriority(entry))
        ]
    }
    
    private def convert(Iterable<String> suffix) {
        suffix.map[idValueConverter.toString(it)].join('.')
    }
    
    private def String getDescription(ILayoutMetaData data) {
        '''«data.name» («data.id»)'''
    }
    
    override protected getCrossrefFilter(CrossReference reference, ContentAssistContext context) {
        val model = context.currentModel
        if (model instanceof ElkEdgeSection) {
            switch reference {
                case grammar.elkEdgeSectionAccess.incomingShapeElkConnectableShapeCrossReference_4_0_0_2_0,
                case grammar.elkSingleEdgeSectionAccess.incomingShapeElkConnectableShapeCrossReference_1_0_0_2_0:
                    return new SectionShapeFilter(model, SectionShapeFilter.INCOMING)
                case grammar.elkEdgeSectionAccess.outgoingShapeElkConnectableShapeCrossReference_4_0_1_2_0,
                case grammar.elkSingleEdgeSectionAccess.outgoingShapeElkConnectableShapeCrossReference_1_0_1_2_0:
                    return new SectionShapeFilter(model, SectionShapeFilter.OUTGOING)
            }
        }
        // Show only cross-references from the same resource
        // (cross-resource references don't work anyway due to the opposite references)
        val resourceURI = model.eResource.URI
        return [ candidate |
            candidate.EObjectURI.trimFragment == resourceURI
        ]
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