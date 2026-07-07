/*******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.grandom.generators;

import java.util.Random;

/**
 * Utility class which generates blind text.
 */
public final class BlindTextGenerator {
	
    /** The blind text the generator is based on. */
    private static final String LOREM_IPSUM = "Lorem ipsum dolor sit amet, consetetur sadipscing elitr,"
            + " sed diam nonumy eirmod tempor invidunt ut labore et dolore magna aliquyam erat,"
            + " sed diam voluptua. At vero eos et accusam et justo duo dolores et ea rebum. "
            + "Stet clita kasd gubergren, no sea takimata sanctus est Lorem ipsum dolor sit amet.";
    /** The generated text will have at least this many characters. */
    private static final int CHARACTER_MIN_LENGTH = 40;
    /** We randomly add up to this many characters onto the minimum text length. */
    private static final int CHARACTER_RANGE = 10;
    
    /**
     * As a utility class, this is not supposed to be instantiated.
     */
    private BlindTextGenerator() {        
    }

    /**
     * Generate a random excerpt of our blind text.
     */
    public static String generate() {
        Random randomGenerator = new Random();
        
        // Choose the number of characters
        int numOfCharacters = randomGenerator.nextInt(CHARACTER_RANGE);
        numOfCharacters += CHARACTER_MIN_LENGTH;
        
        // Choose where to start extracting from our base text
        int offset = randomGenerator.nextInt(LOREM_IPSUM.length() - (CHARACTER_RANGE + CHARACTER_MIN_LENGTH));

        return LOREM_IPSUM.substring(offset, offset + numOfCharacters).trim();
    }

}
