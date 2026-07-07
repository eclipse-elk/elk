/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.data;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.eclipse.elk.core.data.LayoutOptionData.Target;
import org.eclipse.elk.core.data.LayoutOptionData.Visibility;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;

import com.google.common.base.Joiner;
import com.google.common.base.Strings;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

/**
 * Utility class providing basic mechanisms for layout option content assist.
 */
public final class LayoutDataContentAssist {

    private LayoutDataContentAssist() { }

    // SUPPRESS CHECKSTYLE NEXT 30 Javadoc|VisibilityModifier
    /** Container with proposal information. */
    public static final class Proposal<T> {

        /** The proposed string representation to be used. */
        public final String proposal;
        /**
         * A speaking label describing the proposals content. May be {@code null} in which case {@link #proposal} shall
         * be used.
         */
        public final String label;
        /** The underlying object, e.g. an instance of the corresponding {@link LayoutOptionData}. */
        public final T data;

        private Proposal(final String proposal, final String label, final T data) {
            this.proposal = proposal;
            this.label = label;
            this.data = data;
        }

        public static <T> Proposal<T> of(final String proposal, final String label, final T data) {
            return new Proposal<T>(proposal, label, data);
        }

        public static <T> Proposal<T> of(final String proposal, final String label) {
            return new Proposal<T>(proposal, label, null);
        }

        public static <T> Proposal<T> of(final String proposal, final T data) {
            return new Proposal<T>(proposal, null, data);
        }
    }

    /**
     * @return a list of layout option proposals where the layout options' ids match {@code prefix} and are applicable
     *         to {@code element}.
     */
    public static List<Proposal<LayoutOptionData>> getLayoutOptionProposals(final ElkGraphElement element,
            final String prefix) {

        LayoutAlgorithmData responsibleAlgorithm = getAlgorithm(element);
        if (element instanceof ElkNode) {
            Set<Target> targets = ((ElkNode) element).getChildren().isEmpty() ? ImmutableSet.of(Target.NODES)
                    : ImmutableSet.of(Target.NODES, Target.PARENTS);
            return getLayoutOptionProposals(element, responsibleAlgorithm, targets, prefix);
        } else if (element instanceof ElkEdge) {
            return getLayoutOptionProposals(element, responsibleAlgorithm, Target.EDGES, prefix);
        } else if (element instanceof ElkPort) {
            return getLayoutOptionProposals(element, responsibleAlgorithm, Target.PORTS, prefix);
        } else if (element instanceof ElkLabel) {
            return getLayoutOptionProposals(element, responsibleAlgorithm, Target.LABELS, prefix);
        }

        return Collections.emptyList();
    }

    /**
     * @return a list of layout algorithm proposals where the layout algorithms' ids match {@code prefix}.
     */
    public static List<Proposal<LayoutAlgorithmData>> getLayoutAlgorithmProposals(final String prefix) {
        return LayoutMetaDataService.getInstance().getAlgorithmData().stream()
                .map(o -> matchesPrefix(o, LayoutMetaDataService.getInstance()::getAlgorithmDataBySuffix, prefix))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.toList());
    }
    
    /**
     * @return a list of layout option value proposals, where the layout {@code option}'s value matches {@code prefix}.
     */
    public static List<Proposal<Object>> getLayoutOptionValueProposal(final LayoutOptionData option,
            final String prefix) {
        List<Proposal<Object>> proposals = Lists.newArrayList();

        switch (option.getType()) {
        case BOOLEAN:
        case ENUM:
        case ENUMSET:
            String[] choices = option.getChoices();
            for (int i = 0; i < choices.length; ++i) {
                String proposal = choices[i];
                String label = proposal;
                Enum<?> enumVal = option.getEnumValue(i);
                if (ElkGraphUtil.isExperimentalPropertyValue(enumVal)) {
                    label += " - Experimental";
                } else if (ElkGraphUtil.isAdvancedPropertyValue(enumVal)) {
                    label += " - Advanced";
                }
                proposals.add(Proposal.of(proposal, label));
            }
            break;
        case DOUBLE:
        case INT:
            proposals.add(Proposal.of(option.getDefaultDefault().toString(), option.getType().toString()));
            break;
        case OBJECT:
            String proposal = "";
            try {
                proposal = "\"" + option.getOptionClass().newInstance().toString() + "\"";
            } catch (InstantiationException e) {
                // fail silent
            } catch (IllegalAccessException e) {
                // fail silent
            }
            proposals.add(Proposal.of(proposal, option.getType().toString()));
            break;
        default:
            // nothing to add
        }

        return proposals;
    }

    private static List<Proposal<LayoutOptionData>> getLayoutOptionProposals(final ElkGraphElement element,
            final LayoutAlgorithmData algorithmData, final LayoutOptionData.Target targetType, final String prefix) {
        return getLayoutOptionProposals(element, algorithmData,
                targetType != null ? Sets.newHashSet(targetType) : Collections.emptySet(), prefix);
    }

    private static List<Proposal<LayoutOptionData>> getLayoutOptionProposals(final ElkGraphElement element,
            final LayoutAlgorithmData algorithmData, final Set<LayoutOptionData.Target> targetTypes,
            final String prefix) {
        return LayoutMetaDataService.getInstance().getOptionData().stream()
                // Gather all layout options that come into question
                .filter(o -> (o.getTargets() == null || o.getTargets().isEmpty())
                        || !Collections.disjoint(targetTypes, o.getTargets()))
                .filter(o -> algorithmData == null || algorithmData.knowsOption(o) || o.equals(CoreOptions.ALGORITHM))
                .filter(o -> element == null || !element.getProperties().containsKey(o))
                .filter(o -> o.getVisibility() != Visibility.HIDDEN)
                // Match them against the passed prefix
                .map(o -> matchesPrefix(o, LayoutMetaDataService.getInstance()::getOptionDataBySuffix, prefix))
                // Finally collect
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.toList());
    }

    /**
     * See {@link LayoutDataContentAssistTest} for what this method is expected to accept.
     * 
     * What do we want to accept?
     * 
     * <pre>
     * 1) Empty prefix
     * 2) Any substring of the complete id that starts after a dot, e.g. 'org.eclipse.elk.algorithm', 'elk.algorithm'
     * 3) Any substring of the layout option name, e.g. "lgo" for "Layout Algorithm"
     * </pre>
     */
    private static <T extends ILayoutMetaData> Optional<Proposal<T>> matchesPrefix(final T data,
            final Function<String, T> checker, final String prefix) {

        // Case 3)
        final boolean matchesName =
                !Strings.isNullOrEmpty(prefix) && data.getName().toLowerCase().contains(prefix.toLowerCase());

        List<String> idSplit = Arrays.asList(data.getId().split("\\."));
        List<String> prefixSplit = matchesName ? null : Arrays.asList(prefix.split("\\."));

        // For layout options, guarantee that the proposed option always starts with the group
        int start = idSplit.size() - 1;
        if (data instanceof LayoutOptionData) {
            LayoutOptionData layoutOption = (LayoutOptionData) data;
            if (!Strings.isNullOrEmpty(layoutOption.getGroup())) {
                --start;
            }
            long numberOfGroups = layoutOption.getGroup().chars().filter(c -> c == '.').count();
            start -= numberOfGroups;
        }

        for (int i = start; i >= 0; i--) {
            List<String> suffixElements = idSplit.subList(i, idSplit.size());
            String suffix = Joiner.on('.').join(suffixElements);
            // Case 2) (note that Case 1 is included here)
            if (checker.apply(suffix) != null && (matchesName || startsWith(suffixElements, prefixSplit))) {
                return Optional.of(Proposal.of(suffix, data));
            }
        }

        return Optional.empty();
    }

    /**
     * Checks if {@code strings} (if it were joined in some way) at some starts with {@code prefix}. To give an example:
     * {@code [elk, spacing, nodeNode]} starts both with {@code [spacing, nodeN]} and {@code [nodeN]}.
     */
    private static boolean startsWith(final List<String> strings, final List<String> prefix) {
        if (prefix.isEmpty()) {
            return true;
        }
        for (int i = 0; i < strings.size() - prefix.size() + 1; i++) {
            int j = 0;
            boolean matches = true;
            while (j < prefix.size() && matches) {
                matches = strings.get(i + j).startsWith(prefix.get(j));
                j++;
            }
            if (matches) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return the algorithm that is responsible for the layout of the given element. Note that this might be ambiguous:
     *         e.g. a port of a composite node can be handled both by the algorithm that arranges that node and the
     *         algorithm that arranges its container node.
     */
    public static LayoutAlgorithmData getAlgorithm(final ElkGraphElement element) {
        ElkNode relevantNode = getRelevantNode(element);
        if (relevantNode != null) {
            String algorithmId = relevantNode.getProperty(CoreOptions.ALGORITHM);
            if (!Strings.isNullOrEmpty(algorithmId)) {
                return LayoutMetaDataService.getInstance().getAlgorithmDataBySuffix(algorithmId);
            }
        }
        return null;
    }
    
    private static ElkNode getRelevantNode(final ElkGraphElement element) {
        if (element instanceof ElkNode) {
            ElkNode node = (ElkNode) element;
            if (node.getParent() != null) {
                return node.getParent();
            } else {
                return node;
            }
        } else if (element instanceof ElkEdge) {
            ElkEdge edge = (ElkEdge) element;
            return edge.getContainingNode();
        } else if (element instanceof ElkPort) {
            ElkPort port = (ElkPort) element;
            if (port.getParent() != null && port.getParent().getParent() != null) {
                return port.getParent().getParent();
            } else if (port.getParent() != null) {
                return port.getParent();
            }
        } else if (element instanceof ElkLabel) {
            ElkLabel label = (ElkLabel) element;
            ElkGraphElement parent = label.getParent();
            while (parent instanceof ElkLabel) {
                ElkLabel parentLabel = (ElkLabel) parent;
                parent = parentLabel.getParent();
            }
            return getRelevantNode(parent);
        }

        return null;
    }
    
}
