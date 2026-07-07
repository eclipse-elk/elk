/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
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
package org.eclipse.elk.graph.json.text.ide.contentassist

import com.google.inject.Inject
import com.google.inject.Provider
import org.eclipse.elk.core.data.ILayoutMetaData
import org.eclipse.elk.core.data.LayoutDataContentAssist
import org.eclipse.elk.core.data.LayoutOptionData
import org.eclipse.elk.core.data.LayoutOptionData.Type
import org.eclipse.elk.core.options.CoreOptions
import org.eclipse.elk.graph.ElkGraphElement
import org.eclipse.elk.graph.impl.ElkPropertyToValueMapEntryImpl
import org.eclipse.elk.graph.json.text.services.ElkGraphJsonGrammarAccess
import org.eclipse.xtext.Alternatives
import org.eclipse.xtext.Assignment
import org.eclipse.xtext.Keyword
import org.eclipse.xtext.RuleCall
import org.eclipse.xtext.conversion.impl.IDValueConverter
import org.eclipse.xtext.ide.editor.contentassist.ContentAssistContext
import org.eclipse.xtext.ide.editor.contentassist.ContentAssistEntry
import org.eclipse.xtext.ide.editor.contentassist.IIdeContentProposalAcceptor
import org.eclipse.xtext.ide.editor.contentassist.IdeContentProposalProvider
import org.eclipse.xtext.util.Strings

/**
 * Special content assist proposals for the ELK Graph JSON language.
 */
class ElkGraphJsonProposalProvider extends IdeContentProposalProvider {
    
    static val DISABLED_KEYWORDS = #{'}', ']'}

    ElkGraphJsonGrammarAccess grammar
    IDValueConverter idValueConverter

    @Inject
    def void initialize(Provider<IDValueConverter> idValueConverterProvider, ElkGraphJsonGrammarAccess grammarAccess) {
        this.idValueConverter = idValueConverterProvider.get => [
            rule = grammarAccess.IDRule
        ]
        this.grammar = grammarAccess
    }

    override protected filterKeyword(Keyword keyword, ContentAssistContext context) {
        !DISABLED_KEYWORDS.contains(keyword.value) && keyword.value != context.prefix
    }

    override protected _createProposals(Keyword keyword, ContentAssistContext context,
        IIdeContentProposalAcceptor acceptor) {
        if (filterKeyword(keyword, context)) {
            val entry = proposalCreator.createProposal(keyword.value, context, ContentAssistEntry.KIND_KEYWORD, null)
            if (entry !== null) {
                entry.kind = ContentAssistEntry.KIND_KEYWORD
                entry.source = keyword
                acceptor.accept(entry, proposalPriorities.getKeywordPriority(keyword.value, entry))
            }
        }
    }
    
    override protected _createProposals(RuleCall ruleCall, ContentAssistContext context,
        IIdeContentProposalAcceptor acceptor) {
        // By default rule calls are not handled. 
        // The following code ensures that the known keys (e.g. 'x', 'ports') are properly 
        // suggested. Although only the variant without any quotes is suggested. 
        if (ruleCall.rule.name.startsWith("Key")) {
            switch alternatives: ruleCall.rule.alternatives {
                Alternatives:
                    alternatives.elements
                        .filter(Keyword)
                        .filter[!it.value.startsWith('"') && !it.value.startsWith("'")]
                        .forEach[it.createProposals(context, acceptor)]
            }
        }
    }
    
    override protected _createProposals(Assignment assignment, ContentAssistContext context,
        IIdeContentProposalAcceptor acceptor) {
        switch assignment {
            case grammar.propertyAccess.keyAssignment_0:
                completePropertyKey(context, acceptor)
            case grammar.propertyAccess.valueAssignment_2_0,
            case grammar.propertyAccess.valueAssignment_2_1,
            case grammar.propertyAccess.valueAssignment_2_2:
                completePropertyValue(context, acceptor)
            default:
                super._createProposals(assignment, context, acceptor)
        }
    }
    
    protected def void completePropertyKey(ContentAssistContext context, IIdeContentProposalAcceptor acceptor) {
        val model = context.currentModel
        if (model instanceof ElkGraphElement) {
            LayoutDataContentAssist.getLayoutOptionProposals(model, context.prefix.unquoteIfNecessary).forEach [ p |
                val entry = new ContentAssistEntry => [
                    proposal = p.proposal.convertPropertyId
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
    
    /** Besides {@link CoreOptions#ALGORITHM}, there's another option that allows to select a layout algorithm. 
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
    
    protected def proposeAlgorithms(ContentAssistContext context, IIdeContentProposalAcceptor acceptor) {
        LayoutDataContentAssist.getLayoutAlgorithmProposals(context.prefix.unquoteIfNecessary).forEach [ p |
            val entry = new ContentAssistEntry => [
                proposal = p.proposal.convertPropertyId
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
    
    private def typeAwarePropertyValueProposal(LayoutOptionData option, ContentAssistContext context,
        IIdeContentProposalAcceptor acceptor) {
        LayoutDataContentAssist.getLayoutOptionValueProposal(option, context.prefix.unquoteIfNecessary).forEach [ p |
            val proposal = if (#{Type.ENUM, Type.ENUMSET, Type.STRING}.contains(option.type)) {
                    '"' + p.proposal + '"'
                } else {
                    p.proposal
                }
            val entry = proposalCreator.createProposal(proposal, context, ContentAssistEntry.KIND_VALUE) [
                label = p.label
                source = option
            ]
            acceptor.accept(entry, proposalPriorities.getDefaultPriority(entry))
        ]
    }
    
    private def convertPropertyId(String proposal) {
        return '"' + Strings.split(proposal, '.').map[idValueConverter.toString(it)].join('.') + '"'
    }

    private def String getDescription(ILayoutMetaData data) {
        '''«data.name» («data.id»)'''
    }

    private def unquoteIfNecessary(String s) {
        var unquoted = s
        while (unquoted.startsWith('"') || unquoted.startsWith("'")) {
            unquoted = unquoted.substring(1, unquoted.length)
        }
        return unquoted
    }
    
}