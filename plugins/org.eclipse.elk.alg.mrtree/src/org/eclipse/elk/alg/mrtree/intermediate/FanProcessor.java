/*******************************************************************************
 * Copyright (c) 2013, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree.intermediate;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Map.Entry;

import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * A processor that computes the maximal fan out for each node in the given graph. The maximal fan out of
 * a node is the maximal number of descendants it has got in one level.
 * 
 * @author sor
 * @author sgu
 */
public class FanProcessor implements ILayoutProcessor<TGraph> {

    /** this map temporarily contains the maximal fan of each node in the given graph. */
    private Map<String, Integer> gloFanMap = new HashMap<String, Integer>();

    /** this map temporarily contains the number of descendants of each node in the given graph. */
    private Map<String, Integer> gloDescMap = new HashMap<String, Integer>();

    @Override
    public void process(final TGraph tGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Processor compute fanout", 1);

        /** clear map if processor is reused */
        gloFanMap.clear();
        gloDescMap.clear();

        TNode root = null;
        /** find the root of the component */
        Iterator<TNode> it = tGraph.getNodes().iterator();
        while (root == null && it.hasNext()) {
            TNode tNode = it.next();
            if (tNode.getProperty(InternalProperties.ROOT)) {
                root = tNode;
            }
        }

        LinkedList<TNode> rootLevel = new LinkedList<TNode>();
        rootLevel.add(root);

        calculateFan(rootLevel);

        /** set the fan and descendants for all nodes */
        for (TNode tNode : tGraph.getNodes()) {
            String key = tNode.getProperty(InternalProperties.ID);
            int fan = gloFanMap.get(key) != null ? gloFanMap.get(key) : 0;
            tNode.setProperty(InternalProperties.FAN, fan);
            int desc = 1 + (gloDescMap.get(key) != null ? gloDescMap.get(key) : 0);
            tNode.setProperty(InternalProperties.DESCENDANTS, desc);
        }

        progressMonitor.done();
    }

    /**
     * Calculates a fan for the nodes in the list currentLevel including their descendants.
     * 
     * @param currentLevel
     *            the list of TNode for which the fan should be calculated.
     */
    private void calculateFan(final LinkedList<TNode> currentLevel) {
        if (!currentLevel.isEmpty()) {

            /** the children of the current level are the next level */
            LinkedList<TNode> nextLevel = new LinkedList<TNode>();

            /** currentLevel is not empty so stringId will be set */
            String id = null;
            String pId = null;

            /** the size of the which the stringId will be extended for this level */
            int digits = (int) (Math.floor(Math.log10(currentLevel.size())) + 1);

            /**
             * set final stringId for all nodes in this level and set the provisional stringId for
             * their children
             */
            int index = 0;
            for (TNode tNode : currentLevel) {
                /** the final stringId is the stringId of the parent and the extension */
                if (pId != tNode.getProperty(InternalProperties.ID)) {
                    pId = tNode.getProperty(InternalProperties.ID);
                    index = 0;
                }
                if (pId != null) {
                    id = pId + formatRight(index++, digits);
                } else {
                    id = formatRight(index++, digits);
                }
                tNode.setProperty(InternalProperties.ID, id);
                for (TNode tChild : tNode.getChildren()) {
                    nextLevel.add(tChild);
                    /** the provisional stringId is the Id of the parent */
                    tChild.setProperty(InternalProperties.ID, id);
                }
            }

            /** holds the occurences of descendants in this level */
            Map<String, Integer> locFanMap = new HashMap<String, Integer>();

            /** calculated occurences of descendants in this level */
            for (int i = 0; i < id.length() - digits; i++) {
                for (TNode tNode : currentLevel) {
                    String key = tNode.getProperty(InternalProperties.ID).substring(0, i + 1);
                    int blockSize = locFanMap.get(key) != null ? locFanMap.get(key) + 1 : 1;
                    locFanMap.put(key, blockSize);
                }
            }

            /** update the gloFanMap with the values from locFanMap will only increase the values */
            for (Entry<String, Integer> locEntry : locFanMap.entrySet()) {
                Integer gloValue = gloDescMap.get(locEntry.getKey()) != null ? gloDescMap
                        .get(locEntry.getKey()) : 0;
                gloDescMap.put(locEntry.getKey(), locEntry.getValue() + gloValue);
                gloValue = gloFanMap.get(locEntry.getKey());
                if (gloValue == null || gloValue < locEntry.getValue()) {
                    gloFanMap.put(locEntry.getKey(), locEntry.getValue());
                }
            }

            /** calculate the occurrences in the deeper levels and add them to the global gloFanMap */
            calculateFan(nextLevel);
        }
    }
    
    /**
     * Create a string representation of the given value, padded with 0s.
     * 
     * @param value a value
     * @param len the min. number of characters
     * @return a string representation
     */
    public static String formatRight(final int value, final int len) {
        String s = value + "";
        while (s.length() < len) {
            s = "0" + s;
        }
        return s;
    }
    
}
